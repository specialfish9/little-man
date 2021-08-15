package mnkgame.players;

import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;
import mnkgame.*;

/*
 * Minimax + NegaScout (Alpha Beta variation for ID) + Iterative Deepening + k-1 heuristics
 * + euristic evaluation at iterative depth + caching (memoization) via zobrist hashing
 */
public class GalileoGalilei implements MNKPlayer {
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
      key = nextZobrist(i, j);
      return result;
    }

    @Override
    public void unmarkCell() {
      MNKCell last = MC.get(MC.size() - 1);
      key = nextZobrist(last.i, last.j);
      super.unmarkCell();
    }

    // computes the new/old hash (xor works both ways)
    public long nextZobrist(int i, int j) {
      return key ^ zobrist[i * M + j][currentPlayer];
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

  private static final double HALT = Double.MIN_VALUE;
  private static final double MIN = -Double.MAX_VALUE, MAX = Double.MAX_VALUE;
  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, ENEMY_WIN;
  private int M, N, K;
  private long startTime, timeout;
  private Random r;
  private Board board;

  // transposition table hashed using zobrist method
  private HashMap<Long, Pair<Double, Integer>> cache = new HashMap<>();
  private AtomicBoolean zobristReady = new AtomicBoolean(false);
  private static long[][] zobrist;
  private static int maxDepth;
  private static MNKCell bestCell;

  // NOTE: profiling
  private static int visited, cutoff, seriesFound, evaluated, cacheHits, cacheMisses;
  private static boolean clear = true, verbose = true;

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
                    playerName()
                        + "\t: cache ready, in another thread after "
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
    if (cc == null) return null;
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

  private double mem(double result, int depth) {
    if (zobristReady.get()) cache.put(board.zobrist(), new Pair<>(result, depth));
    return result;
  }

  private double unmem(int depth) {
    long hash = board.zobrist();
    if (zobristReady.get() && cache.containsKey(hash)) {
      Pair<Double, Integer> val = cache.get(hash);
      if (val.second == depth) return val.first;
      else return HALT;
    } else return HALT;
  }

  // selection sort for an array of MNKCells and relative evaluation data which
  // places just one minimum/maximum at the beginning of the array.
  // Has therefore a cost of O(n) and is simpler than the quick select
  // algorithm with the median of medians partition function.
  public void selectionSort(Action action, MNKCell[] vec, double[] values, int start) {
    int m = start;
    // find the min in [start,end]
    for (int i = start + 1; i < vec.length; i++)
      if (action == Action.MAXIMIZE ? (values[i] < values[m]) : (values[i] > values[m])) m = i;

    // if we found a new min/max put it in start
    if (m != start) {
      MNKCell tmp = vec[start];
      vec[m] = vec[start];
      vec[start] = tmp;
    }
  }

  private double minimax(Action action, int depth, Double a, Double b) {
    if (board.gameState() != MNKGameState.OPEN) {
      double value = evaluate();
      int actualDepth = Math.max(board.getFreeCells().length - depth, 1);
      return mem(action == Action.MAXIMIZE ? value / actualDepth : value * actualDepth, depth);
    }
    // TODO: check conflicts
    double cached;
    if ((cached = unmem(depth)) != HALT) {
      cacheHits++;
      return cached;
    } else cacheMisses++;

    if (depth == 0) {
      double value = evaluate();
      int actualDepth = Math.max(board.getFreeCells().length - depth, 1);
      return mem(action == Action.MAXIMIZE ? value / actualDepth : value * actualDepth, depth);
    }
    if (shouldHalt()) return HALT;

    visited++;

    MNKCell omc = null; // one move win/loss cell
    if (board.getMarkedCells().length >= 2 * K - 1
        && ((omc = findOneMoveWin(action == Action.MAXIMIZE ? MY_WIN : ENEMY_WIN)) != null
            || (omc = findOneMoveLoss(action == Action.MAXIMIZE ? ENEMY_WIN : MY_WIN)) != null)) {
      seriesFound++; // found once cell in a k-1 series

      board.markCell(omc.i, omc.j);
      // NOTE: must memoize here as oterwhise we'll use the unmarked board hash
      double result = mem(minimax(opposite(action), depth - 1, a, b), depth - 1);
      board.unmarkCell();
      return result;
    }

    MNKCell bc = null;
    MNKCell[] cells = board.getFreeCells();
    double[] ratings = new double[cells.length];
    for (int i = 0; i < cells.length; i++) {
      long hash = board.nextZobrist(cells[i].i, cells[i].j);
      ratings[i] = cache.containsKey(hash) ? cache.get(hash).first : 0;
    }

    int i = 0;
    while (i < cells.length) {
      selectionSort(action, cells, ratings, i);
      MNKCell c = cells[i];
      board.markCell(c.i, c.j);

      double score;
      if (i == 0) score = -minimax(opposite(action), depth - 1, -b, -a);
      else {
        score = -minimax(opposite(action), depth - 1, -a - 1, -a);
        if (score == HALT) {
          board.unmarkCell();
          break;
        }

        if (score > a && score < b) score = -minimax(opposite(action), depth - 1, -b, -score);
      }
      board.unmarkCell();

      if (score == HALT) break;

      if (score > a) { // a = max(a, score)
        a = score;
        bc = c;
      }
      if (a >= b) {
        cutoff += board.getFreeCells().length - 1 - i;
        break;
      }

      i++; // go onto the next move
    }
    if (depth == maxDepth) bestCell = bc;
    if (a == MAX || a == MIN) return HALT;
    else return mem(a, depth);
  }

  public MNKCell iterativeDeepening() {
    int len = board.getFreeCells().length;
    double result = MIN;
    MNKCell bc = null;

    MNKCell omc = null; // one move win/loss cell
    if (board.getMarkedCells().length >= (K - 1) * 2
        && ((omc = findOneMoveWin(MY_WIN)) != null || (omc = findOneMoveLoss(ENEMY_WIN)) != null)) {
      seriesFound++; // found once cell in a k-1 series
      return omc;
    }

    // as long as we can deepen the tree for this move
    // {{{ deepening
    maxDepth = 1;
    while (!shouldHalt()) {
      double latest = minimax(Action.MINIMIZE, maxDepth, MIN, MAX);
      if (latest == HALT) break;

      result = latest; // always keep the last result as it's the most accurate
      bc = bestCell;

      // we can breack when we have a tuple of the type (x, x, false) which means
      // the minimax result was completely obtained by non-heuristics values
      // we can also break when we find a state which will grant us certain win
      if (result >= winCutoff) break;

      maxDepth++;
    }
    // }}}

    // TODO: remove in production
    if (verbose)
      System.out.println(
          playerName() + "\t: finished at depth " + maxDepth + " with value: " + result);

    return bc;
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    startTime = System.currentTimeMillis();
    visited = cutoff = seriesFound = evaluated = cacheHits = cacheMisses = 0;

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
      System.out.println(playerName() + "\t: cut off " + cutoff + " branches");
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
    return "Galileo Galilei";
  }
}
