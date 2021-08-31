package mnkgame.players;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import java.util.Stack;
import java.util.concurrent.atomic.AtomicBoolean;
import mnkgame.*;

public class GalileoGalilei implements MNKPlayer {
  // {{{ tuple
  private static class Pair<A, B> {
    public A first;
    public B second;

    public Pair(A first, B second) {
      this.first = first;
      this.second = second;
    }

    @Override
    public String toString() {
      return "(" + first + "," + second + ")";
    }
  }

  private static class Tuple<A, B, C> {
    public A first;
    public B second;
    public C third;

    public Tuple(A first, B second, C third) {
      this.first = first;
      this.second = second;
      this.third = third;
    }

    @Override
    public String toString() {
      return "(" + first + "," + second + "," + third + ")";
    }
  }
  // }}}

  // {{{ board
  private static class Board extends MNKBoard {
    private int minMN;
    private final MNKCellState me;
    private long key = 0;
    private int queueP1 = 0, queueP2 = 0, queueFree = 0;
    private Queue<MNKCellState> queue = new LinkedList<>();
    private Stack<Integer> previousValues = new Stack<>();
    private int value = 0;

    public Board(int M, int N, int K, int minMN, MNKCellState me) {
      super(M, N, K);
      this.minMN = minMN;
      this.me = me;
    }

    // returns the depth of this board relative to the root of the game where the
    // board was empty. It is equal to the length of the marked cells.
    public int depth() {
      return MC.size();
    }

    public MNKGameState markCell(int i, int j, boolean updateValue) {
      // mind the order of the calls
      key = nextZobrist(i, j);
      if (updateValue) {
        previousValues.push(value);
        value -= eval(i, j);
      }
      MNKGameState result = super.markCell(i, j);
      if (updateValue) value += eval(i, j);
      return result;
    }

    @Override
    public MNKGameState markCell(int i, int j) {
      return markCell(i, j, true);
    }

    public void unmarkCell(boolean revertValue) {
      // mind the order of the calls
      MNKCell last = MC.getLast();
      super.unmarkCell();
      key = nextZobrist(last.i, last.j);
      if (revertValue) value = previousValues.pop();
    }

    @Override
    public void unmarkCell() {
      unmarkCell(true);
    }

    // computes the hash for a new mark (xor works both ways, but pay attention to
    // the
    // currentPlayer)
    public long nextZobrist(int i, int j) {
      return key ^ zobrist[i * minMN + j][(currentPlayer + 1) % 2];
    }

    public int value() {
      return value;
    }

    public long zobrist() {
      return key;
    }

    private int eval(int i, int j) {
      int value = 0;

      // column
      queueClear();
      // for (int ii = Math.max(i - K, 0); ii <= Math.min(i + K, M - 1); ii++)
      for (int ii = 0; ii < M; ii++) value += pushCell(B[ii][j]);

      // row
      queueClear();
      for (int jj = 0; jj < N; jj++) value += pushCell(B[i][jj]);

      // diagonal
      int ku = Math.min(i, j),
          kl = Math.min(M - 1 - i, N - 1 - j),
          ii = i - ku,
          jj = j - ku,
          iim = i + kl,
          jjm = j + kl;
      queueClear();
      for (; ii <= iim && jj <= jjm; ii++, jj++) value += pushCell(B[ii][jj]);

      // counter diagonal
      ii = i - ku;
      jj = j + ku;
      iim = i + kl;
      jjm = j - kl;
      queueClear();
      for (; ii <= iim && jj <= jjm; ii++, jj--) value += pushCell(B[ii][jj]);

      return value;
    }

    private void queueClear() {
      queueP1 = queueP2 = queueFree = 0;
      queue.clear();
    }

    private void popCell() {
      MNKCellState state = queue.poll();
      if (state == MNKCellState.FREE) queueFree--;
      else if (state == MNKCellState.P1) queueP1--;
      else if (state == MNKCellState.P2) queueP2--;
    }

    private int pushCell(final MNKCellState state) {
      if (queue.size() >= K) // useless >
      popCell();

      if (state == MNKCellState.FREE) queueFree++;
      else if (state == MNKCellState.P1) queueP1++;
      else if (state == MNKCellState.P2) queueP2++;
      queue.add(state);
      int sign = me == MNKCellState.P1 ? 1 : -1;
      if (queueP1 + queueFree == K) return sign * (seriesValue(queueFree) + (queueP1*queueP1));
      else if (queueP2 + queueFree == K) return -sign * (seriesValue(queueFree) + (queueP2*queueP2));
      else return 0;
    }

    private int seriesValue(final int free) {
      if (K > 4 && free == 3) return 1000; // 1k
      if (K > 3 && free == 2) return 100000; // 100k
      if (K > 2 && free == 1) return 10000000; // 10M

      return 0;
    }

    @Override
    public String toString() {
      String str = "";
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++)
          str +=
              " " + (B[i][j] == MNKCellState.P1 ? 'x' : (B[i][j] == MNKCellState.P2 ? 'o' : '-'));
        str += '\n';
      }
      return str;
    }
  }
  // }}}

  // {{{ transposition cleanup
  private class CleanupRunnable implements Runnable {
    private final long endTime;
    private final MNKCell move;

    public CleanupRunnable(long e, MNKCell m) {
      endTime = e;
      move = m;
    }

    private boolean shouldHalt() {
      return System.currentTimeMillis() >= endTime;
    }

    private void remove(MNKCell root) {
      // prevent from removing any boards that have the chosen move
      // marked as these boards could come in handy later
      if (root == move || shouldHalt()) return;

      cacheBoard.markCell(root.i, root.j, false);
      if (cache.containsKey(board.zobrist())) {
        cache.remove(board.zobrist());
        if (cacheBoard.gameState() == MNKGameState.OPEN)
          for (MNKCell node : cacheBoard.getFreeCells()) remove(node);
      }
      cacheBoard.unmarkCell(false);
    }

    public void run() {
      for (MNKCell c : cacheBoard.getFreeCells()) {
        if (Thread.currentThread().isInterrupted() || shouldHalt()) break;

        // TODO: recursive
        remove(c);
      }
      System.out.println(playerName() + "\t: cleanup ended");
    }
  }

  private Thread cleanupThread = null;

  private void cleanup(long e, MNKCell m) {
    stopCleanup();
    cleanupThread = new Thread(new CleanupRunnable(e, m));
    cleanupThread.start();
  }

  private void stopCleanup() {
    if (cleanupThread != null && !cleanupThread.isAlive()) cleanupThread.interrupt();
    try {
      cleanupThread.join();
    } catch (Exception e) {
    }
    cleanupThread = null;
  }
  // }}}

  private static final int INFTY = 1000000000; // 1B
  private static final int HALT = -INFTY * 2; // -2B
  private static final double SAFETY_THRESHOLD = 0.9;

  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, ENEMY_WIN;
  private int M, N, K, minMN;
  private long startTime, timeout;
  private Random r;
  private Board board, cacheBoard;

  // transposition table hashed using zobrist method
  private HashMap<Long, int[]> cache = new HashMap<>();
  private AtomicBoolean zobristReady = new AtomicBoolean(false);
  private static long[][] zobrist;

  // NOTE: profiling
  private static int minimaxed,
      failed,
      failedDepth,
      cutoff,
      seriesFound,
      evaluated,
      cacheHits,
      cacheMisses;
  private static boolean clear = false,
      verbose = true,
      superVerbose = false,
      winCutoffEnabled = true;

  public void initPlayer(int M, int N, int K, boolean first, int timeoutInSecs) {
    this.M = M;
    this.N = N;
    this.K = K;
    this.minMN = Math.min(M, N);
    timeout = timeoutInSecs * 1000;
    startTime = System.currentTimeMillis();
    MY_WIN = first ? MNKGameState.WINP1 : MNKGameState.WINP2;
    ENEMY_WIN = first ? MNKGameState.WINP2 : MNKGameState.WINP1;
    ME = first ? MNKCellState.P1 : MNKCellState.P2;
    ENEMY = first ? MNKCellState.P2 : MNKCellState.P1;

    r = new Random(startTime);
    board = new Board(M, N, K, minMN, ME);
    cacheBoard = new Board(M, N, K, minMN, ME);
    cache.clear();

    // continue filling the table for the zobrist hashing function in another
    // thread as to avoid failing initialization
    zobristReady.set(false);
    zobrist = new long[M * N][2];
    int i;
    for (i = 0; i < zobrist.length; i++) {
      if (i % 10 == 0 && shouldHalt()) // check every 10 iterations
      break;

      zobrist[i][0] = r.nextLong();
      zobrist[i][1] = r.nextLong();
    }

    // it did not finish in time, continue in a separate thread
    // NOTE: we can safely modify the zobrist array across threads as it's only
    // read from after the
    if (i != zobrist.length) {
      final int j = i;
      new Thread(
              () -> {
                for (int k = j; k < zobrist.length; k++) {
                  zobrist[k][0] = r.nextLong();
                  zobrist[k][1] = r.nextLong();
                }
                zobristReady.set(true);
                // TODO: remove
                System.out.println(
                    playerName()
                        + "\t: cache ready, in another thread after "
                        + (System.currentTimeMillis() - startTime)
                        + "ms");
              })
          .start();
    } else {
      zobristReady.set(true);
      // TODO: remove
      System.out.println(playerName() + "\t: cache ready, no thread needed");
    }

    // TODO: remove in production
    // clear screen
    if (clear) {
      System.out.print("\033[H\033[2J");
      System.out.flush();
    }
  }

  private boolean shouldHalt() {
    // TODO: tweak values
    return (System.currentTimeMillis() - startTime)
        >= timeout * SAFETY_THRESHOLD; // livin' on the edge
  }

  // finds the first cell needed to copmlete a K-1 streak in any possible direction
  private MNKCell findOneMoveWin(final MNKGameState winState) {
    for (MNKCell c : board.getFreeCells()) {
      MNKGameState result = board.markCell(c.i, c.j);
      board.unmarkCell();
      if (result == winState) return c;
    }
    return null;
  }

  private MNKCell pickRandomNonClosingCell(MNKCell previous) {
    for (MNKCell c : board.getFreeCells()) {
      // avoiding the shouldHalt check here, kinda superflous
      MNKGameState result = board.markCell(c.i, c.j);
      board.unmarkCell();
      if (result == MNKGameState.OPEN
          && (previous == null || previous.i != c.i || previous.j != c.j)) return c;
    }
    return null;
  }

  // finds the first cell the enemy needs to copmlete a K-1 streak in any possible direction
  private MNKCell findOneMoveLoss(final MNKGameState lossState) {
    MNKCell randomCell = null;
    if (board.getFreeCells().length == 1 || (randomCell = pickRandomNonClosingCell(null)) == null)
      return null; // cannot check for enemy's next move when it doesn't exist

    board.markCell(randomCell.i, randomCell.j);
    MNKCell c = findOneMoveWin(lossState);
    board.unmarkCell(); // remove the marked randomCell
    if (c != null) return c;

    // test the randomCell we selected at first. It may be a one-move loss cell
    // get a new random cell different from the previous and call it cc
    MNKCell cc = pickRandomNonClosingCell(randomCell);
    if (board.markCell(cc.i, cc.j)
        != MNKGameState.OPEN) { // randomCell puts us in a draw, ignore that
      board.unmarkCell();
      return null;
    }
    MNKGameState result =
        board.markCell(randomCell.i, randomCell.j); // let the enemy take the random ane
    board.unmarkCell();
    board.unmarkCell();
    return result == lossState ? randomCell : null;
  }

  // evaluate a state (either heuristically or in a deterministic way) regardless
  // of its depth. The depth is also taken into account
  private int evaluate() {
    // check for draws first, most lickely
    MNKGameState state = board.gameState();
    if (state == MNKGameState.DRAW) return 0;
    else if (state == MY_WIN) return (INFTY - 1) / board.depth();
    else if (state == ENEMY_WIN) return (-INFTY + 1) / board.depth();
    else {
      evaluated++;
      return board.value() / board.depth();
    }
  }

  // selection sort for an array of MNKCells and relative evaluation data which
  // places just one minimum/maximum at the beginning of the array.
  // Has therefore a cost of O(n) and is simpler than the quick select
  // algorithm with the median of medians partition function.
  public void selectionSort(MNKCell[] vec, int[] values, int start, int color) {
    int m = start;
    // find the max/min in [start,end]
    for (int i = start + 1; i < vec.length; i++)
      if (color > 0 ? values[i] > values[m] : values[i] < values[m]) m = i;

    // if we found a new min/max put it in start
    if (m != start) {
      MNKCell tmp = vec[m];
      vec[m] = vec[start];
      vec[start] = tmp;

      int tmp1 = values[m];
      values[m] = values[start];
      values[start] = tmp1;
    }
  }

  private Pair<MNKCell[], int[]> getMoves(int searchDepth) {
    MNKCell[] cells = board.getFreeCells();
    int[] ratings = new int[cells.length];
    for (int i = 0; i < cells.length; i++) {
      long hash = board.nextZobrist(cells[i].i, cells[i].j);
      if (cache.containsKey(hash)) {
        int values[] = cache.get(hash);
        ratings[i] = values[1] >= searchDepth - 1 ? values[0] : 0;
      } else ratings[i] = 0;

      // TODO: remove in production
      if (ratings[i] != 0) cacheHits++;
      else cacheMisses++;
    }
    return new Pair<>(cells, ratings);
  }

  // {{{ Principal Variation Search for subtrees
  private int pvs(int color, int searchDepth, int alpha, int beta) {
    int[] entry = {-color * INFTY, searchDepth};
    if (cache.containsKey(board.zobrist())
        && (entry = cache.get(board.zobrist()))[1] >= searchDepth) {
      cacheHits++;
      return entry[0];
    } else cacheMisses++;

    int result;
    MNKCell omc = null;
    if (searchDepth <= 0 || board.gameState() != MNKGameState.OPEN) result = color * evaluate();
    else if (shouldHalt()) result = HALT;
    else if (board.depth() >= 2 * K - 1
        && // only check for one-win-moves
        // if there have been placed
        // enough cells to make one happen
        ((omc = findOneMoveWin(color > 0 ? MY_WIN : ENEMY_WIN)) != null
            || (omc = findOneMoveLoss(color > 0 ? ENEMY_WIN : MY_WIN)) != null)) {
      board.markCell(omc.i, omc.j);
      result = -pvs(-color, searchDepth - 1, -beta, -alpha);
      board.unmarkCell();
    } else {
      if (superVerbose) {
        for (int i = 0; i < board.depth(); i++) System.out.print("--");
        System.out.println("opened " + color + " ( " + alpha + " , " + beta + " )");
      }
      Pair<MNKCell[], int[]> moves = getMoves(searchDepth);
      int bestMove = 0;
      selectionSort(moves.first, moves.second, 0, color);
      board.markCell(moves.first[bestMove].i, moves.first[bestMove].j);
      result = -pvs(-color, searchDepth - 1, -beta, -alpha);
      if (superVerbose) {
        for (int i = 0; i < board.depth(); i++) System.out.print("--");
        System.out.println(" best: " + result);
      }
      board.unmarkCell();
      if (result == HALT) return HALT;
      if (result > alpha) {
        if (result >= beta) return result;
        alpha = result;
      }

      int i = 1, len = board.getFreeCells().length;
      while (i < len) {
        selectionSort(moves.first, moves.second, i, color);
        board.markCell(moves.first[i].i, moves.first[i].j);
        int score = -pvs(-color, searchDepth - 1, -alpha - 1, -alpha); // null window search
        if (score > alpha && score < beta && score != HALT) {
          // research with window [alfa;beta]
          score = -pvs(-color, searchDepth - 1, -beta, -alpha);
          if (score > alpha) alpha = score;
        }
        if (superVerbose) {
          for (int j = 0; j < board.depth(); j++) System.out.print("--");
          System.out.println(" score: " + score);
        }
        board.unmarkCell();
        if (score == HALT) break;
        if (score > result) {
          result = score;
          bestMove = i;
          if (score >= beta) break;
        }
        i++;
      }
      if (superVerbose) {
        for (int j = 0; j < board.depth(); j++) System.out.print("--");
        System.out.println("closed " + result + " ( " + alpha + " , " + beta + " )");
      }
    }
    if (result == HALT) return HALT;

    entry[0] = result;
    cache.put(board.zobrist(), entry);
    return result;
  }
  // }}}

  // {{{ Principal Variation Search on root
  private Pair<Integer, MNKCell> pvsRoot(int searchDepth, int alpha, int beta) {
    int bestValue = -INFTY;
    MNKCell bestCell = null;
    if (board.depth() >= 2 * K - 1
        && // only check for one-win-moves
        // if there have been placed
        // enough cells to make one happen
        ((bestCell = findOneMoveWin(MY_WIN)) != null
            || (bestCell = findOneMoveLoss(ENEMY_WIN)) != null)) {
      board.markCell(bestCell.i, bestCell.j);
      bestValue = -pvs(-1, searchDepth - 1, -beta, -alpha);
      board.unmarkCell();
    } else {
      int bestMove = 0;
      Pair<MNKCell[], int[]> moves = getMoves(searchDepth);
      selectionSort(moves.first, moves.second, 0, 1);
      board.markCell(moves.first[bestMove].i, moves.first[bestMove].j);
      bestValue = -pvs(-1, searchDepth - 1, -beta, -alpha);
      board.unmarkCell();
      if (bestValue == HALT) return new Pair<>(HALT, null);
      if (bestValue > alpha) {
        if (bestValue >= beta) return new Pair<>(bestValue, moves.first[bestMove]);
        alpha = bestValue;
      }

      int i = 1, len = board.getFreeCells().length;
      while (i < len) {
        selectionSort(moves.first, moves.second, i, 1);
        board.markCell(moves.first[i].i, moves.first[i].j);
        int score = -pvs(-1, searchDepth - 1, -alpha - 1, -alpha); // alphaBeta or zwSearch
        if (score > alpha && score < beta && score != HALT) {
          // research with window [alfa;beta]
          score = -pvs(-1, searchDepth - 1, -beta, -alpha);
          if (score > alpha) alpha = score;
        }
        board.unmarkCell();
        if (score == HALT) return new Pair<>(HALT, null);
        if (score > bestValue) {
          bestValue = score;
          bestMove = i;
          if (score >= beta) break;
        }
        i++;
      }
      bestCell = moves.first[bestMove];
    }
    return new Pair<>(bestValue, bestCell);
  }
  // }}}

  // {{{ iterative deepening
  public Pair<Integer, MNKCell> iterativeDeepening() {
    int len = board.getFreeCells().length;
    Pair<Integer, MNKCell> value = null;

    int maxDepth = 1;
    while (!shouldHalt() && maxDepth <= len) {
      // the alpha and beta values are given by the highest and lowest possible
      // achievable values with our evaluation function. Because of how it's coded
      // (as it takes into account the absolute depth of the board) we can compute
      // a much tighter maximum value for alpha/beta and use this to achieve higher
      // cutoffs
      int max = INFTY / (board.depth() + maxDepth);
      Pair<Integer, MNKCell> latest = pvsRoot(maxDepth, -max, max);
      if (latest.first == HALT || latest.first == -HALT) break;

      value = latest;
      // stop the search if we found a certain win
      if (winCutoffEnabled && value.first >= max) break;

      // TODO: remove in production
      if (verbose)
        System.out.println("minimax went to depth " + maxDepth + " with value: " + value);
      maxDepth++;
    }

    // TODO: remove in production
    if (maxDepth - 1 > board.getFreeCells().length) {
      throw new Error(
          "FATAL: iterativeDeepening exceeded allowable depth with a depth of: " + maxDepth);
    }

    return value;
  }
  // }}}

  // {{{ selectCell
  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    startTime = System.currentTimeMillis();
    minimaxed =
        failed = failedDepth = cutoff = seriesFound = evaluated = cacheHits = cacheMisses = 0;
    stopCleanup();

    if (MC.length > 0) {
      board.markCell(
          MC[MC.length - 1].i, MC[MC.length - 1].j); // keep track of the opponent's marks
      if (MC.length > 1)
        cacheBoard.markCell(
            MC[MC.length - 2].i, MC[MC.length - 2].j); // keep track of my previous move
      cacheBoard.markCell(
          MC[MC.length - 1].i, MC[MC.length - 1].j); // keep track of the opponent's marks
    }

    try {
      Pair<Integer, MNKCell> result = iterativeDeepening();
      System.out.println(
          playerName()
              + "\t: visited "
              + minimaxed
              + " nodes, ended with result: ("
              + result.first
              + ", "
              + result.second
              + ")");
      System.out.println(
          playerName()
              + "\t: cut off "
              + cutoff
              + " branches ("
              + (cutoff / Math.max(minimaxed, 1) * 100)
              + "%). failed pvs "
              + failed
              + " times with an average depth of "
              + (failedDepth * 1d) / failed);
      System.out.println(
          playerName()
              + "\t: found a total of "
              + seriesFound
              + " free cells in series (up to k-1)");
      System.out.println(playerName() + "\t: heuristically evaluated " + evaluated + " boards");
      int perc = (int) ((1d * cacheHits / (cacheMisses + cacheHits)) * 100);
      System.out.println(
          playerName()
              + "\t: cached "
              + cache.size()
              + " elements, hit: "
              + cacheHits
              + ", misses: "
              + cacheMisses
              + ". rate: "
              + perc
              + "%");
      System.out.println(
          playerName() + "\t: took " + (System.currentTimeMillis() - startTime) + "ms");

      // TODO: remove in prod
      if (FC.length != board.getFreeCells().length) {
        System.out.println("FATAL: minimax didn't clean the board");
        return FC[0];
      }

      // TODO: remove in production
      if (result.second == null) {
        System.out.println("FATAL: no chosen cell");
        result.second = FC[new Random().nextInt(FC.length)];
      }

      if (board.markCell(result.second.i, result.second.j) == MNKGameState.OPEN)
        cleanup(System.currentTimeMillis() + (long) (timeout * SAFETY_THRESHOLD), result.second);
      return result.second;
    } catch (Exception e) {
      e.printStackTrace();
      MNKCell c = FC[new Random().nextInt(FC.length)];
      cleanup(System.currentTimeMillis() + (long) (timeout * SAFETY_THRESHOLD), c);
      board.markCell(c.i, c.j);
      return c;
    }
  }
  // }}}

  public String playerName() {
    return "Galileo Galilei";
  }
}
