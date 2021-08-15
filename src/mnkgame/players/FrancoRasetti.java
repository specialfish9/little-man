package mnkgame.players;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;
import mnkgame.*;

/*
 * Minimax + Alpha Beta + Iterative Deepening + k-1 heuristics
 * + euristic evaluation at iterative depth + caching (memoization)
 */
public class FrancoRasetti implements MNKPlayer {
  // {{{ tuple
  private static class Pair<A extends Comparable<A>, B> {
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

  private static class PairComparator<A extends Comparable<A>, B>
      implements Comparator<Pair<A, B>> {
    @Override
    public int compare(Pair<A, B> c1, Pair<A, B> c2) {
      return c1.first.compareTo(c2.first);
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
    private long key = 0;

    public Board(int M, int N, int K) {
      super(M, N, K);
    }

    @Override
    public MNKGameState markCell(int i, int j) {
      MNKGameState result = super.markCell(i, j);
      updateZobrist(i * M + j, currentPlayer);
      return result;
    }

    @Override
    public void unmarkCell() {
      MNKCell last = MC.get(MC.size() - 1);
      updateZobrist(last.i * M + last.j, currentPlayer);
      super.unmarkCell();
    }

    // computes the new/old hash (xor works both ways)
    private void updateZobrist(int index, int player) {
      key ^= zobrist[index][player];
    }

    public long zobrist() {
      return key;
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

  // {{{ chances
  // TODO: check per entrambi i giocatori uno come per Series (?)
  private static class Chances {
    private static Board board;
    private static int M, N, K;
    private static double klog;
    private static MNKCellState player;

    private static int currStreak, currStreakFree;

    private static double checkCell(int i, int j) {
      MNKCellState c = board.cellState(i, j);
      if (c == player || c == MNKCellState.FREE) {
        currStreak++;
        currStreakFree++;
      } else currStreak = currStreakFree = 0;

      double value = 0;
      if (currStreak >= K) {
        // we know for a fact that currStreakFree > 0, oterwhise this method
        // won't be called as the evaluation can be done in a deterministic way
        value = Math.log(currStreak / currStreakFree) / klog;
        // should be 1 when currStreakFree == 1 which is the best value and
        // the greatest chance to win (only need to mark 1 cell) whereas
        // should get closer to 0 as currStreakFree increases. It's effectively
        // log_k(currStreak/currStreakFree)
        currStreak--;
        currStreakFree--;
      }
      return value;
    }

    private static void reset() {
      currStreak = currStreakFree = 0;
    }

    // Returns the number of possible series the player could streak to win the game
    public static double winningChances(final Board b, final MNKCellState p) {
      board = b;
      player = p;
      M = board.M;
      N = board.N;
      K = board.K;
      klog = Math.log(K);

      double chance = 0;

      // check all columns
      for (int i = 0; i < M && M - i + currStreak >= K; i++) {
        reset();
        for (int j = 0; j < N && N - j + currStreak >= K; j++) chance += checkCell(i, j);
      }

      // check all rows
      for (int j = 0; j < N && N - j + currStreak >= K; j++) {
        reset();
        for (int i = 0; i < M && M - i + currStreak >= K; i++) chance += checkCell(i, j);
      }

      // iterate over all diagonals
      int nDiagonals = (Math.min(N, M) - K) * 2 + 1;
      for (int x = 0; x < nDiagonals; x++) {
        reset();
        int i = 0, j = 0;
        if (x != 0 && x % 2 == 0) i = x / 2;
        else if (x != 0) j = x;

        // TODO: don't check useless cells like in prev loops
        while (i < M && j < N) {
          checkCell(i, j);
          j++;
          i++;
        }
      }

      // iterate over all counter diagonals
      for (int x = 0; x < nDiagonals; x++) {
        reset();
        int i = 0, j = N - 1;
        if (x != 0 && x % 2 == 0) i = x / 2;
        else if (x != 0) j = N - 1 - x;

        while (i < M && j >= 0) {
          checkCell(i, j);
          j--;
          i++;
        }
      }

      return chance;
    }
  }
  // }}}

  // {{{ search result
  private static class SearchResult {
    public double value;
    public int depth;
    public boolean heuristic;

    public SearchResult(double value, int depth, boolean heuristic) {
      this.value = value;
      this.depth = depth;
      this.heuristic = heuristic;
    }

    @Override
    public String toString() {
      return "(value=" + value + ",depth=" + depth + ",heuristic=" + heuristic + ")";
    }
  }
  // }}}

  private static final double HALT = Double.MIN_VALUE;
  private static final double MIN = -Double.MAX_VALUE, MAX = Double.MAX_VALUE;
  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, ENEMY_WIN;
  private int M, N, K;
  private long startTime, timeout;
  private Random r;
  private Board board;

  // transposition table hashed using zobrist method
  private HashMap<Long, SearchResult> cache = new HashMap<>();
  private AtomicBoolean zobristReady = new AtomicBoolean(false);
  private static long[][] zobrist;

  // NOTE: profiling
  private static int visited, seriesFound, evaluated, cacheHits, cacheMisses;
  private static boolean clear = false, verbose = false;

  public void initPlayer(int M, int N, int K, boolean first, int timeoutInSecs) {
    this.M = M;
    this.N = N;
    this.K = K;
    timeout = timeoutInSecs;
    startTime = System.currentTimeMillis();
    MY_WIN = first ? MNKGameState.WINP1 : MNKGameState.WINP2;
    ENEMY_WIN = first ? MNKGameState.WINP2 : MNKGameState.WINP1;
    ME = first ? MNKCellState.P1 : MNKCellState.P2;
    ENEMY = first ? MNKCellState.P2 : MNKCellState.P1;

    r = new Random(startTime);
    board = new Board(M, N, K);
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
      Thread t =
          new Thread(
              () -> {
                for (int k = j; k < zobrist.length; k++) {
                  zobrist[k][0] = r.nextLong();
                  zobrist[k][1] = r.nextLong();
                }
                zobristReady.set(true);
                // TODO: remove
                System.out.println(
                    playerName() + "\t: cache ready, in another thread after "
                        + (System.currentTimeMillis() - startTime) / 1000d
                        + "s");
              });
      t.start();
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

  enum Action {
    MINIMIZE,
    MAXIMIZE
  }

  private static Action opposite(final Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  private boolean shouldHalt() {
    // TODO: tweak values
    return (System.currentTimeMillis() - startTime) / 1000.0
        >= timeout * 0.90; // livin' on the edge
  }

  private boolean shouldHalt(long endTime) {
    return System.currentTimeMillis() >= endTime;
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
    if (board.markCell(cc.i, cc.j) != MNKGameState.OPEN) {
      // randomCell puts us in a draw, ignore that
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
  // of its depth. The depth will be taken into account later to allow for caching
  private static final double winCutoff =
      2d; // represents a certain win, and therefore we can cutoff search

  private double evaluate() {
    // check for draws first, most lickely
    if (board.gameState() == MNKGameState.DRAW) return 0;
    else if (board.gameState() == MY_WIN) return 2 * M * N; // 2 times max depth
    else if (board.gameState() == ENEMY_WIN) return -2;
    else {
      evaluated++;
      // keep the heuristic evaluation between 1 and -1
      double res =
          Math.min(
              Math.max(
                  (int) (Chances.winningChances(board, ME) - Chances.winningChances(board, ENEMY)),
                  -1),
              1);
      return res;
    }
  }

  private SearchResult mem(SearchResult result) {
    if (zobristReady.get()) cache.put(board.zobrist(), result);
    return result;
  }

  private SearchResult unmem() {
    if (zobristReady.get())
      return cache.containsKey(board.zobrist()) ? cache.get(board.zobrist()) : null;
    else return null;
  }

  private SearchResult minimax(Action action, int depth, long endTime, Double a, Double b) {
    if (board.gameState() != MNKGameState.OPEN) {
      double value = evaluate();
      int actualDepth = Math.max(board.getFreeCells().length - depth, 1);
      return mem(
          new SearchResult(
              action == Action.MAXIMIZE ? value / actualDepth : value * actualDepth, depth, false));
    }
    // TODO: check conflicts
    SearchResult cached;
    if ((cached = unmem()) != null && cached.depth == depth) {
      cacheHits++;
      return cached;
    } else cacheMisses++;

    if (depth == 0) {
      double value = evaluate();
      int actualDepth = Math.max(board.getFreeCells().length - depth, 1);
      return mem(
          new SearchResult(
              action == Action.MAXIMIZE ? value / actualDepth : value * actualDepth, depth, true));
    }
    if (shouldHalt(endTime)) return new SearchResult(HALT, depth, true);

    visited++;

    MNKCell omc = null; // one move win/loss cell
    if (board.getMarkedCells().length >= (K - 1) * 2
        && // only check for one-win-moves
        // if there have been placed
        // enough cells to make one happen
        ((omc = findOneMoveWin(action == Action.MAXIMIZE ? MY_WIN : ENEMY_WIN)) != null
            || (omc = findOneMoveLoss(action == Action.MAXIMIZE ? ENEMY_WIN : MY_WIN)) != null)) {
      seriesFound++; // found once cell in a k-1 series

      board.markCell(omc.i, omc.j);
      // NOTE: must memoize here as oterwhise we'll use the unmarked board hash
      SearchResult result = mem(minimax(opposite(action), depth - 1, endTime, a, b));
      board.unmarkCell();
      return result;
    }

    double best = action == Action.MAXIMIZE ? MIN : MAX;
    boolean isBestHeuristic = false;
    // TODO: (?) shuffle(series.noseries);
    for (MNKCell c : board.getFreeCells()) {
      board.markCell(c.i, c.j);
      SearchResult result = minimax(opposite(action), depth - 1, endTime, a, b);
      board.unmarkCell();

      if (result.value == HALT) break;

      if ((action == Action.MAXIMIZE && result.value > best)
          || (action == Action.MINIMIZE && result.value < best)
          || (result.value == best && isBestHeuristic && !result.heuristic)) {
        best = result.value;
        isBestHeuristic = result.heuristic;
      }
      if (action == Action.MAXIMIZE && best > a) a = best;
      else if (action == Action.MINIMIZE && best < b) b = best;

      if (b <= a) break;
    }
    if (best == MAX || best == MIN) return new SearchResult(HALT, depth, true);
    else return mem(new SearchResult(best, depth, isBestHeuristic));
    // TODO: careful about the third value. Should it be like this? when we find
    // a configuration we like that is not heuristic we stop. But maybe others have
    // been evaluated wrongly with heurisitc, recieved a lower score, and were ignored
    // They could still be valuable. We wanna investigate that.
  }

  public MNKCell iterativeDeepening() {
    int len = board.getFreeCells().length;
    double best = MIN;
    boolean isBestHeuristic = false;
    MNKCell bestCell = null;

    MNKCell omc = null; // one move win/loss cell
    if (board.getMarkedCells().length >= (K - 1) * 2
        && ((omc = findOneMoveWin(MY_WIN)) != null || (omc = findOneMoveLoss(ENEMY_WIN)) != null)) {
      seriesFound++; // found once cell in a k-1 series
      return omc;
    }

    for (MNKCell c : board.getFreeCells()) {
      long endTime = System.currentTimeMillis() + (timeout * 900) / len;
      board.markCell(c.i, c.j);
      int depth = 1;
      SearchResult result = new SearchResult(0, 1, false);

      // as long as we can deepen the tree for this move
      // {{{ deepening
      while (!shouldHalt(endTime)) {
        SearchResult latest = minimax(Action.MINIMIZE, depth, endTime, MIN, MAX);
        if (latest.value == HALT) break;

        result = latest; // always keep the last result as it's the most accurate

        // we can breack when we have a tuple of the type (x, x, false) which means
        // the minimax result was completely obtained by non-heuristics values
        // we can also break when we find a state which will grant us certain win
        if (!result.heuristic || result.value >= winCutoff) break;

        depth++;
      }
      // }}}

      // TODO: remove in production
      if (depth - 1 > board.getFreeCells().length) {
        throw new Error(
            "FATAL: iterativeDeepening exceeded allowable depth with a depth of: " + depth);
      }

      // TODO: remove in production
      if (verbose)
        System.out.println(
            "minimax for cell "
                + board.getMarkedCells()[board.getMarkedCells().length - 1]
                + " finished at depth "
                + depth
                + " with value: "
                + result.value
                + ", heuristic: "
                + result.heuristic);
      board.unmarkCell();

      // stop the search if we found a winning node
      if (result.value >= winCutoff) return c;

      if (result.value > best || (result.value == best && isBestHeuristic && !result.heuristic)) {
        best = result.value;
        isBestHeuristic = result.heuristic;
        bestCell = c;
      }
    }
    return bestCell;
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    startTime = System.currentTimeMillis();
    visited = seriesFound = evaluated = cacheHits = cacheMisses = 0;

    if (MC.length > 0)
      board.markCell(
          MC[MC.length - 1].i, MC[MC.length - 1].j); // keep track of the opponent's marks

    try {
      MNKCell result = iterativeDeepening();
      System.out.println(
          playerName()
              + "\t: visited "
              + visited
              + " nodes, ended with result: (x,x,"
              + result
              + ")");
      System.out.println(
          playerName()
              + "\t: found a total of "
              + seriesFound
              + " free cells in series (up to k-3)");
      System.out.println(playerName() + "\t: heuristically evaluated " + evaluated + " boards");
      int perc = (int) (((double) cacheHits / (cacheMisses + cacheHits)) * 100);
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
          playerName() + "\t: took " + (System.currentTimeMillis() - startTime) / 1000d + "s");

      // TODO: remove in prod
      if (FC.length != board.getFreeCells().length) {
        System.out.println("FATAL: minimax didn't clean the board");
        return FC[0];
      }

      board.markCell(result.i, result.j);
      return result;
    } catch (Exception e) {
      e.printStackTrace();
      return FC[new Random().nextInt(FC.length)];
    }
  }

  public String playerName() {
    return "Franco Rasetti";
  }
}
