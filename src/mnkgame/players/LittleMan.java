package mnkgame.players;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.Stack;
import java.util.concurrent.atomic.AtomicBoolean;
import mnkgame.*;

// NOTE: the code uses vim's marker folding. Use `za` to toggle folds.

public class LittleBoy implements MNKPlayer {
  private static final int INFTY = 1000000000; // 1 Billion
  // Out of bounds value used to denote a timeout in the search
  private static final int HALT = -INFTY * 2; // -2 Billion
  // NOTE: tweak as needed to prevent exceeding time limits.
  private static final double SAFETY_THRESHOLD = 0.95;

  // Transposition table keys are hashed using the Zobrist technique
  // Cache entry structure: [marked, lastCell, searchDepth, type, value]
  // Marked and lastCell are used to avoid collisions. Type can be one of:
  // EXACT_VALUE => value is the exact evaluation of the board
  // UPPER_BOUND => value is the upper bound
  // LOWER_BOUND => value is the lower bound
  private static final int EXACT_VALUE = 0, UPPER_BOUND = 1, LOWER_BOUND = -1;
  private HashMap<Long, int[]> cache = new HashMap<>();
  // Zobrist key fragments generated randomly during init
  private long[][] zobrist;
  private AtomicBoolean isZobristReady = new AtomicBoolean(false);
  private MNKGameState MY_WIN, ENEMY_WIN;
  private int M, N, K, minMN;
  private Random r;
  private Board board;

  private long startTime, timeout;
  // Return true when we're close to running out of time (as defined by the
  // SAFETY_THRESHOLD constant).
  private boolean shouldHalt() {
    return (System.currentTimeMillis() - startTime) >= timeout;
  }

  public String playerName() {
    return "Little Boy";
  }

  // {{{ board

  // An extension of the provided MNKBoard to account for incremental hash
  // generation and board evaluation. All methods have been implemented in such
  // a way to maintain the same asymptotic cost of the original implementation.
  private class Board extends MNKBoard {
    // Zobrist hash value
    private long key = 0;

    // Heuristic value and previous values (for fast unmakes)
    private Stack<Integer> previousValues = new Stack<>();
    private int value = 0;

    private final boolean first;
    private final int minMN;

    public Board(int M, int N, int K, int minMN, boolean first) {
      super(M, N, K);
      this.minMN = minMN;
      this.first = first;
    }

    // Returns the depth of this board relative to the root.
    public int marked() {
      return MC.size();
    }

    @Override
    public MNKGameState markCell(int i, int j) {
      return markCell(i, j, true);
    }

    // Marks the given cell at (i, j) via callin the super methods, but also
    // computes the new Zobrist key and the new statc value of the match.
    // Cost: O(max{M, N}) because of `markCell` and `eval`
    public MNKGameState markCell(int i, int j, boolean updateInternals) {
      // mind the order of the calls
      if (updateInternals) {
        key = nextZobrist(i, j);
        previousValues.push(value);
        value -= eval(i, j);
      }
      MNKGameState result = super.markCell(i, j);
      if (updateInternals) value += eval(i, j);
      return result;
    }

    @Override
    public void unmarkCell() {
      unmarkCell(true);
    }

    // Ummarks the last move and restores the previous hash and evaluation values.
    // Cost: O(1)
    public void unmarkCell(boolean updateInternals) {
      // mind the order of the calls
      MNKCell last = MC.getLast();
      super.unmarkCell();
      if (updateInternals) {
        key = nextZobrist(last.i, last.j);
        value = previousValues.pop();
      }
    }

    // Computes the hash for a new mark. We can use this same method to both do
    // and undo the hash value, thanks to the properties of the XOR gate.
    // Care must be taken _when_ calling this as it relies on the currentPlayer
    // being correctly set to the one who has done/is undoing the move.
    public long nextZobrist(int i, int j) {
      return key ^ zobrist[i * minMN + j][(currentPlayer + 1) % 2];
    }

    // Returns the Zobrist hash value for the current board
    public long zobrist() {
      return key;
    }

    private int queueP1 = 0, queueP2 = 0, queueFree = 0;

    // Assigns aribrary high values to k-3, k-2, and k-1 series to make them
    // stand out from others and help the search move towards creating longer
    // streaks
    private int seriesValue(final int free) {
      if (K > 4 && free == 3) return 1000; // 1k
      if (K > 3 && free == 2) return 100000; // 100k
      if (K > 2 && free == 1) return 10000000; // 10M

      return 0;
    }
    
    // Returns wheter the given (i,j) position is within the bounds of the
    // current board
    private boolean isInBounds(int i, int j) {
      return i >= 0 && i < M && j >= 0 && j < N;
    }
    
    private int pushCell(final int i, final int j, final int dI, final int dJ) {
      if (queueFree + queueP1 + queueP2 >= K) { // useless >
        MNKCellState s = B[i-dI*(K-1)][j-dJ*(K-1)];
        if (s == MNKCellState.FREE) queueFree--;
        else if (s == MNKCellState.P1) queueP1--;
        else if (s == MNKCellState.P2) queueP2--;
      }

      int freeFactor = 1;
      if (B[i][j] == MNKCellState.FREE) queueFree++;
      else if (B[i][j] == MNKCellState.P1) queueP1++;
      else if (B[i][j] == MNKCellState.P2) queueP2++;

      if(queueP1 + queueFree == K || queueP2 + queueFree == K) {
        // Check if the cell after the series is free
        if(isInBounds(i+dI, j+dJ) && B[i+dI][j+dJ] == MNKCellState.FREE)
          freeFactor *= 2;
        // Check if the cell before the series is free
        // NOTE: isInBounds is pedantic here
        if(isInBounds(i-dI*K, j-dJ*K) && B[i-dI*K][j-dJ*K] == MNKCellState.FREE)
          freeFactor *= 2;

        int sign = freeFactor * (first ? 1 : -1);
        if (queueP1 + queueFree == K)
          return sign * (seriesValue(queueFree) + (queueP1 * queueP1));
        else
          return -sign * (seriesValue(queueFree) + (queueP2 * queueP2));
      } else return 0;
    }

    private void queueClear() {
      queueP1 = queueP2 = queueFree = 0;
    }

    private int eval(int i, int j) {
      int value = 0;

      // column
      queueClear();
      // for (int ii = Math.max(i - K, 0); ii <= Math.min(i + K, M - 1); ii++)
      for (int ii = 0; ii < M; ii++) value += pushCell(ii, j, 1, 0);

      // row
      queueClear();
      for (int jj = 0; jj < N; jj++) value += pushCell(i, jj, 0, 1);

      // diagonal
      int ku = Math.min(i, j),
          kl = Math.min(M - i, N - j),
          ii = i - ku,
          jj = j - ku,
          iim = i + kl,
          jjm = j + kl;
      queueClear();
      for (; ii < iim && jj < jjm; ii++, jj++) value += pushCell(ii, jj, 1, 1);

      // counter diagonal
      ii = i - ku;
      jj = j + ku;
      iim = i + kl;
      jjm = j - kl;
      for (; ii < iim && jj < jjm; ii++, jj--) value += pushCell(ii, jj, 1, -1);

      return value;
    }

    public int value() {
      return value;
    }
  }
  // }}}

  // {{{ transposition cleanup

  // The reference to the currently running cleanup thread
  private Thread cleanupThread = null;

  // Runnable which cleans up any cached boards with less marked cells than the
  // current amount + 1 (the one played by the enemy). These boards are
  // effectively useless as MiniMax won't reach them again and they only
  // increase the chance of having cache conflicts.
  private class CleanupRunnable implements Runnable {
    private final long endTime;
    private final int marked;

    public CleanupRunnable(long endTime, int marked) {
      this.endTime = endTime;
      this.marked = marked;
    }

    private boolean shouldHalt() {
      return System.currentTimeMillis() >= endTime;
    }

    // Iterate over all cached entries in our time limit and remove any board
    // which is deemed no longer needed.
    public void run() {
      Iterator<Map.Entry<Long, int[]>> iter = cache.entrySet().iterator();
      while (iter.hasNext()) {
        if (Thread.currentThread().isInterrupted() || shouldHalt()) break;

        Map.Entry<Long, int[]> e = iter.next();
        if (e.getValue()[0] <= marked + 1) iter.remove();
      }
    }
  }

  // Starts the cache cleanup in another thread for the given time limit and amount
  // of currently marked cells.
  private void cleanup(long endTime, int marked) {
    stopCleanup();
    cleanupThread = new Thread(new CleanupRunnable(endTime, marked));
    cleanupThread.start();
  }

  // Stops the cleanup forcefully by sending an interrupt to the thread. We have
  // no interest in stopping gracefully and want this to complete as fast as possible.
  private void stopCleanup() {
    if (cleanupThread != null && cleanupThread.isAlive()) cleanupThread.interrupt();
    try {
      cleanupThread.join();
    } catch (Exception e) {
    }
    cleanupThread = null;
  }
  // }}}

  // {{{ init
  public void initPlayer(int M, int N, int K, boolean first, int timeoutInSecs) {
    this.M = M;
    this.N = N;
    this.K = K;
    minMN = Math.min(M, N);

    MY_WIN = first ? MNKGameState.WINP1 : MNKGameState.WINP2;
    ENEMY_WIN = first ? MNKGameState.WINP2 : MNKGameState.WINP1;

    r = new Random(startTime);
    board = new Board(M, N, K, minMN, first);
    timeout = (long) (timeoutInSecs * 1000 * SAFETY_THRESHOLD);
    startTime = System.currentTimeMillis();
    stopCleanup();
    cache.clear();

    // If necessary, continue filling the table for the zobrist hashing function
    // in another thread to avoid failing initialization.
    isZobristReady.set(false);
    zobrist = new long[M * N][2];
    int i;
    for (i = 0; i < zobrist.length; i++) {
      if (i % 10 == 0 && shouldHalt()) // check every 10 iterations
      break;

      zobrist[i][0] = r.nextLong();
      zobrist[i][1] = r.nextLong();
    }

    // If the previous loop didn't manage to fill the table in time we continue
    // the process in a separate thread. We can safely modify the zobrist array
    // across threads as it's only read from after the `isZobristReady` field
    // gets set to true
    if (i != zobrist.length-1) {
      final int j = i;
      new Thread(
              () -> {
                for (int k = j; k < zobrist.length; k++) {
                  zobrist[k][0] = r.nextLong();
                  zobrist[k][1] = r.nextLong();
                }
                isZobristReady.set(true);
              })
          .start();
    } else isZobristReady.set(true);
  }
  // }}}

  // {{{ one-cell threats

  // Looks at all cells and finds the one which completes a k-1 series to achieve
  // the given `winState`. Cost: O(M*N)
  private MNKCell findOneMoveWin(final MNKGameState winState) {
    for (MNKCell c : board.getFreeCells()) {
      MNKGameState result = board.markCell(c.i, c.j, false);
      board.unmarkCell(false);
      if (result == winState) return c;
    }
    return null;
  }

  // Tries any available cell and returns it if it doesn't change the
  // immediate outcome of the game. Cost: O(M*N), in practice it'll most often
  // do less operations.
  private MNKCell pickRandomNonClosingCell(MNKCell previous) {
    for (MNKCell c : board.getFreeCells()) {
      MNKGameState result = board.markCell(c.i, c.j, false);
      board.unmarkCell(false);
      if (result == MNKGameState.OPEN
          && (previous == null || previous.i != c.i || previous.j != c.j)) return c;
    }
    return null;
  }

  // Finds a possible cell with which the enemy can complete a k-1 series and win
  // Cost: O(M*N)
  private MNKCell findOneMoveLoss(final MNKGameState lossState) {
    MNKCell randomCell = null;
    // Can't check two moves ahead when there aren't enough moves
    if (board.getFreeCells().length <= 2 || (randomCell = pickRandomNonClosingCell(null)) == null)
      return null;

    board.markCell(randomCell.i, randomCell.j, false);
    MNKCell c = findOneMoveWin(lossState);
    board.unmarkCell(false);
    if (c != null) return c;

    // Test the randomCell we selected at first. It may be a one-move loss cell.
    // Get a new random cell, different from the previous, and try with that
    MNKCell cc = pickRandomNonClosingCell(randomCell);
    if (cc == null) return null;
    // TODO: remove, should be granted from the condition above
    // if (board.markCell(cc.i, cc.j, false) != MNKGameState.OPEN) {
    //   // The newrandomCell puts us in a draw, ignore that
    //   board.unmarkCell(false);
    //   return null;
    // }
    board.markCell(cc.i, cc.j, false);
    // Look at the result of the enemy marking the initial random cell
    MNKGameState result = board.markCell(randomCell.i, randomCell.j, false); 
    board.unmarkCell(false);
    board.unmarkCell(false);
    return result == lossState ? randomCell : null;
  }
  // }}}

  // {{{ moves ordering

  // Evaluate a state, either heuristically or in a deterministic way.
  // Cost: O(1)
  private int evaluate() {
    MNKGameState state = board.gameState();
    if (state == MNKGameState.DRAW) return 0;
    else if (state == MY_WIN) return INFTY / board.marked();
    else if (state == ENEMY_WIN) return -INFTY / board.marked();
    else return Math.min(Math.max(board.value(), -(INFTY / 10)), INFTY / 10) / board.marked();
  }

  // Swaps vec[a] with vec[b]. Cost: O(1)
  private <T> void swap(T[] vec, int a, int b) {
    if (a == b) return;

    T tmp = vec[a];
    vec[a] = vec[b];
    vec[b] = tmp;
  }

  // Swaps vec[a] with vec[b]. Cost: \Tehta(1)
  // Copy of the above with native type
  private void swap(int[] vec, int a, int b) {
    int tmp = vec[a];
    vec[a] = vec[b];
    vec[b] = tmp;
  }

  // Single-step Selection Sort for an array of MNKCells and relative evaluation
  // scores. It only applies one step of the Selection Sort algorithm, positioning
  // the minimum/maximum (based on the color) of the vector in the `start` position.
  // Cost: O(end-start)
  public void selectionSort(MNKCell[] vec, int[] values, int start, int end, int color) {
    int m = start;
    // find the max/min in [start,end]
    for (int i = start + 1; i < end; i++)
      if (color > 0 ? values[i] > values[m] : values[i] < values[m]) m = i;

    // swap vec[m] with vec[start] if we found a new max/min
    if(m != start) {
      swap(vec, start, m);
      swap(values, start, m);
    }
  }

  // Selects a random item in the vector [start,end] and swaps it with the
  // item found at the `start` position. Cost: O(1)
  public void randomSelection(MNKCell[] vec, int start, int end) {
    int i = start + r.nextInt(end - start);
    if (i != start) swap(vec, start, i);
  }

  private int rateMoves(MNKCell[] cells, int[] ratings, int searchDepth) {
    int j = cells.length, i = 0; // limits for the [a,b] set containing all not-yet-looked-at cells
    while (i < j) {
      int entry[] =
          cacheEntry(
              board.nextZobrist(cells[i].i, cells[i].j),
              board.marked() + 1,
              cells[i].i * minMN + cells[i].j,
              searchDepth - 1);
      if (entry[3] != 2) {
        ratings[i] = entry[4];
        i++;
      } else {
        swap(cells, i, j - 1);
        ratings[j - 1] = 0;
        j--;
      }
    }
    return j;
  }
  // }}}

  // {{{ Principal Variation Search for subtrees
  private int[] cacheEntry(int searchDepth) {
    MNKCell[] c = board.getMarkedCells();
    return cacheEntry(
        board.zobrist(),
        board.marked(),
        c[c.length - 1].i * minMN + c[c.length - 1].j,
        searchDepth);
  }

  // returns a cache entry for the current board. If the current board is already
  // cached the entry contains the proper data, otherwhise the entry fields 2,3 are
  // dummy. A non-cached board can be identified by entry[3] == 2
  private int[] cacheEntry(long hash, int marked, int lastCell, int searchDepth) {
    if (isZobristReady.get() && cache.containsKey(hash)) {
      int[] cached = cache.get(hash);
      // Make sure the board has the same number of marked symbols and the last
      // cell marked matches. This is done to avoid false positives in the cache
      if (cached[0] == marked
          && cached[1] == lastCell
          && cached[2] >= searchDepth
          // useless
          && cached[3] != 2) return cached;
    }

    return new int[] {marked, lastCell, searchDepth, 2, -INFTY};
  }

  private int pvs(int color, int depth, int alpha, int beta) {
    MNKCell c;
    int prevAlpha = alpha, value = -INFTY;

    // Look up for the current board inside the transposition table. As described
    // in `cacheEntry(4)` entry[2] == 2 notes a cache miss.
    int[] entry = cacheEntry(depth);
    if (entry[3] != 2) {
      if (entry[3] == EXACT_VALUE) return entry[4];
      else if (entry[3] == LOWER_BOUND) alpha = Math.max(alpha, entry[4]);
      else beta = Math.max(beta, entry[4]);

      if (alpha >= beta) return entry[4];
    }

    if (shouldHalt()) return HALT;
    else if (depth <= 0 || board.gameState() != MNKGameState.OPEN) return color * evaluate();
    else if (board.marked() >= 2 * K - 1
        && ((c = findOneMoveWin(color > 0 ? MY_WIN : ENEMY_WIN)) != null
            || (c = findOneMoveLoss(color > 0 ? ENEMY_WIN : MY_WIN)) != null)) {
      board.markCell(c.i, c.j);
      value = -pvs(-color, depth - 1, -beta, -alpha);
      board.unmarkCell();
    } else {
      MNKCell[] moves = board.getFreeCells();
      int[] ratings = new int[moves.length];
      int sortUpTo = rateMoves(moves, ratings, depth);
      for (int i = 0; i < moves.length; i++) {
        if (i < sortUpTo) selectionSort(moves, ratings, i, sortUpTo, color);
        else randomSelection(moves, i, moves.length);

        // NOTE: alpha is only updated when we make a proper full window search
        // to avoid wrong bounds.
        board.markCell(moves[i].i, moves[i].j);
        int score;
        if (i == 0) {
          score = -pvs(-color, depth - 1, -beta, -alpha);
          alpha = Math.max(alpha, score);
        } else {
          // Try first a null window search on non-PV nodes with bounds [-alpha-1, -alpha]
          score = -pvs(-color, depth - 1, -alpha - 1, -alpha);

          // If the search failed inside the [alpha, beta] bounds the result may
          // be meaningful so we need to do a proper search
          if (score > alpha && score < beta && value != HALT) {
            score = -pvs(-color, depth - 1, -beta, -alpha);
            alpha = Math.max(alpha, score);
          }
        }
        board.unmarkCell();

        if (score > value || score == HALT || score == -HALT) value = score;
        if (value >= beta || value == HALT || value == -HALT) break;
      }
    }
    if (value == HALT) return HALT;

    entry[2] = depth;
    entry[4] = value;
    if (value <= prevAlpha) entry[3] = UPPER_BOUND;
    else if (value >= beta) entry[3] = LOWER_BOUND;
    else entry[3] = EXACT_VALUE;

    cache.put(board.zobrist(), entry);
    return value;
  }
  // }}}

  // {{{ Principal Variation Search on root

  private int rootValue = -INFTY;

  private MNKCell pvsRoot(int depth, int alpha, int beta) {
    MNKCell cell = null;
    int value = -INFTY;

    if (shouldHalt()) return null;
    else if (board.marked() >= 2 * K - 1
        && ((cell = findOneMoveWin(MY_WIN)) != null
            || (cell = findOneMoveLoss(ENEMY_WIN)) != null)) {
      board.markCell(cell.i, cell.j);
      value = -pvs(-1, depth - 1, -beta, -alpha);
      board.unmarkCell();
    } else {
      MNKCell[] moves = board.getFreeCells();
      int[] ratings = new int[moves.length];
      int sortUpTo = rateMoves(moves, ratings, depth);
      for (int i = 0; i < moves.length; i++) {
        if (i < sortUpTo) selectionSort(moves, ratings, i, sortUpTo, 1);
        else randomSelection(moves, i, moves.length);

        // NOTE: alpha is only updated when we make a proper full window search
        // to avoid wrong bounds.
        board.markCell(moves[i].i, moves[i].j);
        int score;
        if (i == 0) {
          score = -pvs(-1, depth - 1, -beta, -alpha);
          alpha = Math.max(alpha, score);
        } else {
          // Try first a null window search on non-PV nodes with bounds [-alpha-1, -alpha]
          score = -pvs(-1, depth - 1, -alpha - 1, -alpha);

          // If the search failed inside the [alpha, beta] bounds the result may
          // be meaningful so we need to do a proper search
          if (score > alpha && score < beta && value != HALT) {
            score = -pvs(-1, depth - 1, -beta, -alpha);
            alpha = Math.max(alpha, score);
          }
        }
        board.unmarkCell();
        if (score == HALT || score == -HALT) return null;

        if (score > value) {
          value = score;
          cell = moves[i];
        }
        if (value >= beta) break;
      }
    }

    rootValue = value;
    return cell;
  }
  // }}}

  // {{{ iterative deepening

  // Iterative Deppening calls the pvsRoot depth-first search algorithm with an
  // ever-increasing maximum depth, to search the tree as deep as we can and
  // provide insights about move ordering for future searches.
  // We stop the search if we either timeout or we find a certain win.
  public MNKCell iterativeDeepening() {
    int len = board.getFreeCells().length;
    MNKCell value = null;

    int maxDepth = 1;
    while (!shouldHalt() && maxDepth <= len) {
      // the alpha and beta values are given by the highest and lowest possible
      // achievable values with our evaluation function. Because of how it's coded
      // (as it takes into account the absolute depth of the board) we can compute
      // a much tighter maximum value for alpha/beta and use this to achieve higher
      // cutoffs
      int max = INFTY / Math.min(board.marked() + maxDepth, 2 * K - 1);
      MNKCell latest = pvsRoot(maxDepth, -max, max);
    System.out.println(
        playerName()
            + "\t: at depth "
            + maxDepth
            + " got cell: "
            + latest
            + " with value: "
            + rootValue + " bounds [" + -max + "," + max + "]");

      if (latest == null) break;

      value = latest;

      maxDepth++;
    }

    return value;
  }
  // }}}

  // {{{ selectCell
  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    startTime = System.currentTimeMillis();
    stopCleanup();

    // Keep track of the opponent's marked cells
    if (MC.length > 0)
      board.markCell(MC[MC.length - 1].i, MC[MC.length - 1].j); 

    MNKCell result = iterativeDeepening();
    // to avoid catastrophic failures in case anything breaks
    if (result == null) result = FC[new Random().nextInt(FC.length)];

    if (board.markCell(result.i, result.j) == MNKGameState.OPEN && board.marked() < M * N - 3)
      cleanup(System.currentTimeMillis() + (long) (timeout * SAFETY_THRESHOLD), board.marked());
    return result;
  }
  // }}}
}

// vim: ts=2 sw=2 fdm=marker
