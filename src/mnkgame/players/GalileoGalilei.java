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
    private int minMN;
    private final MNKCellState me;
    private long key = 0;
    private int queueP1 = 0, queueP2 = 0, queueFree = 0;
    private Queue<MNKCellState> queue = new LinkedList<>();
    private Stack<Double> previousValues = new Stack<>();
    private double value = 0;
    private double V[][];

    public Board(int M, int N, int K, int minMN, MNKCellState me) {
      super(M, N, K);
      this.minMN = minMN;
      this.me = me;
      this.V = new double[M][N];
      for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) this.V[i][j] = 0;
    }

    @Override
    public MNKGameState markCell(int i, int j) {
      // mind the order of the calls
      key = nextZobrist(i, j);
      MNKGameState result = super.markCell(i, j);
      double newValue = eval(i, j);
      value += newValue - V[i][j];
      V[i][j] = newValue;
      previousValues.push(value);
      return result;
    }

    @Override
    public void unmarkCell() {
      // mind the order of the calls
      MNKCell last = MC.getLast();
      super.unmarkCell();
      key = nextZobrist(last.i, last.j);
      value = previousValues.pop();
      V[last.i][last.j] = eval(last.i, last.j);
    }

    // computes the hash for a new mark (xor works both ways, but pay attention to
    // the
    // currentPlayer)
    public long nextZobrist(int i, int j) {
      return key ^ zobrist[i * minMN + j][(currentPlayer + 1) % 2];
    }

    // computes the hash for the removal of the last mark (xor works both ways, but
    // we gotta change
    // the currentPlayer)
    public long prevZobrist(int i, int j) {
      return key ^ zobrist[i * minMN + j][currentPlayer];
    }

    public double value() {
      return value;
    }

    public long zobrist() {
      return key;
    }

    private double eval(int i, int j) {
      double value = 0;

      // column
      queueClear();
      for (int ii = Math.max(i - K, 0); ii < Math.min(i + K, M - 1); ii++)
        value += weightCell(B[ii][j]);

      // row
      queueClear();
      for (int jj = Math.max(j - K, 0); jj < Math.min(j + K, N - 1); jj++)
        value += weightCell(B[i][jj]);

      // diagonal
      int ku = Math.min(K, Math.min(i, j)),
          kl = Math.min(K, Math.min(M - 1 - i, N - 1 - j)),
          ii = i - ku,
          jj = j - ku,
          iim = i + kl,
          jjm = j + kl;
      queueClear();
      for (; ii < iim && jj < jjm; ii++, jj++) value += weightCell(B[ii][jj]);

      // TODO: counter diagonal

      return value * (1 / M * N);
    }

    private double weightCell(MNKCellState c) {
      if (queue.size() >= K) // useless >
      queuePop();

      return queuePush(c);
    }

    private void queueClear() {
      queueP1 = queueP2 = queueFree = 0;
      queue.clear();
    }

    private void queuePop() {
      MNKCellState state = queue.poll();
      if (state == MNKCellState.FREE) queueFree--;
      else if (state == MNKCellState.P1) queueP1--;
      else if (state == MNKCellState.P2) queueP2--;
    }

    private double queuePush(final MNKCellState state) {
      if (state == MNKCellState.FREE) queueFree++;
      else if (state == MNKCellState.P1) queueP1++;
      else if (state == MNKCellState.P2) queueP2++;
      queue.add(state);
      if (queue.size() == K) {
        if (queueP1 + queueFree == K) return (me == MNKCellState.P1 ? 1 : -1) * (queueP1 / K);
        else if (queueP2 + queueFree == K) return (me == MNKCellState.P2 ? 1 : -1) * (queueP2 / K);
        else return 0; // TODO: value empty streaks more than filled useless ones
      }
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

  // {{{ chances
  // TODO: check per entrambi i giocatori uno come per Series (?)
  private static class Chances {
    private static Board board;
    private static int M, N, K;
    private static MNKCellState player;

    private static int currStreak, currStreakFree;

    private static double checkCell(int i, int j) {
      MNKCellState c = board.cellState(i, j);
      if (c == player || c == MNKCellState.FREE) {
        currStreak++;
        if (c == MNKCellState.FREE) currStreakFree++;
      } else reset();

      double value = 0;
      if (currStreak >= K && currStreakFree < K) {
        value = 1;
        // TODO: redo (?)
        // we know for a fact that currStreakFree > 0, oterwhise this method
        // won't be called as the evaluation can be done in a deterministic way
        // value = Math.log(currStreak / currStreakFree) / klog;
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

      double chance = 0;

      // check all rows
      for (int i = 0; i < M; i++) {
        reset();
        for (int j = 0; j < N; j++) {
          if (N - j + currStreak < K) break;

          chance += checkCell(i, j);
        }
      }

      // check all columns
      for (int j = 0; j < N; j++) {
        reset();
        for (int i = 0; i < M; i++) {
          if (M - i + currStreak < K) break;

          chance += checkCell(i, j);
        }
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
          chance += checkCell(i, j);
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
          chance += checkCell(i, j);
          j--;
          i++;
        }
      }

      return chance;
    }
  }
  // }}}

  private static final double HALT = Double.MIN_VALUE;
  private static final double MIN = -100, MAX = 2;
  private static final double winCutoff =
      MAX; // represents a certain win, and therefore we can cutoff search

  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, ENEMY_WIN;
  private int M, N, K, minMN, maxDepth = 1;
  private long startTime, timeout;
  private Random r;
  private Board board;

  // transposition table hashed using zobrist method
  private HashMap<Long, double[]> cache = new HashMap<>();
  private AtomicBoolean zobristReady = new AtomicBoolean(false);
  private static long[][] zobrist;

  // NOTE: profiling
  private static int minimaxed, cutoff, seriesFound, evaluated, cacheHits, cacheMisses;
  private static boolean clear = false, verbose = false;

  public void initPlayer(int M, int N, int K, boolean first, int timeoutInSecs) {
    this.M = M;
    this.N = N;
    this.K = K;
    this.minMN = Math.min(M, N);
    timeout = timeoutInSecs;
    startTime = System.currentTimeMillis();
    MY_WIN = first ? MNKGameState.WINP1 : MNKGameState.WINP2;
    ENEMY_WIN = first ? MNKGameState.WINP2 : MNKGameState.WINP1;
    ME = first ? MNKCellState.P1 : MNKCellState.P2;
    ENEMY = first ? MNKCellState.P2 : MNKCellState.P1;

    r = new Random(startTime);
    board = new Board(M, N, K, minMN, ME);
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

  private static Action opposite(final Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  private boolean shouldHalt() {
    // TODO: tweak values
    return (System.currentTimeMillis() - startTime) / 1000.0
        >= timeout * 0.85; // livin' on the edge
  }

  private boolean shouldHalt(long endTime) {
    return System.currentTimeMillis() >= endTime;
  }

  // finds the first cell needed to copmlete a K-1 streak in any possible
  // direction
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

  // finds the first cell the enemy needs to copmlete a K-1 streak in any possible
  // direction
  private MNKCell findOneMoveLoss(final MNKGameState lossState) {
    MNKCell randomCell = null;
    if (board.getFreeCells().length > 2 || (randomCell = pickRandomNonClosingCell(null)) == null)
      return null; // cannot check for enemy's next move when it doesn't exist

    board.markCell(randomCell.i, randomCell.j);
    MNKCell c = findOneMoveWin(lossState);
    board.unmarkCell(); // remove the marked randomCell
    if (c != null) return c;

    // test the randomCell we selected at first. It may be a one-move loss cell
    // get a new random cell different from the previous and call it cc
    MNKCell cc = pickRandomNonClosingCell(randomCell);
    if (cc == null) return null; // cannot find a non-closing cell
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
  // TODO: take depth into account. idea:
  // - certain wins/losses = Math.max(1 / depth, 0.5), Math.max(-0.25*depth, -1);
  // - heuristic values = 0 < win < 0.5, -0.5 < loss < 0
  private double evaluate(int depth) {
    // check for draws first, most lickely
    if (board.gameState() == MNKGameState.DRAW) return 0;
    else if (board.gameState() == MY_WIN) return .5d + 1 / depth;
    else if (board.gameState() == ENEMY_WIN) return -1 - (-.25 * depth);
    else {
      evaluated++;
      // keep the heuristic evaluation between 1 and -1
      return board.value();
    }
  }

  // selection sort for an array of MNKCells and relative evaluation data which
  // places just one minimum/maximum at the beginning of the array.
  // Has therefore a cost of O(n) and is simpler than the quick select
  // algorithm with the median of medians partition function.
  public void selectionSort(MNKCell[] vec, double[] values, int start, int color) {
    int m = start;
    // find the max/min in [start,end]
    for (int i = start + 1; i < vec.length; i++)
      if (color > 0 ? (values[i] > values[m]) : (values[i] < values[m])) m = i;

    // if we found a new min/max put it in start
    if (m != start) {
      MNKCell tmp = vec[m];
      vec[m] = vec[start];
      vec[start] = tmp;

      double tmp1 = values[m];
      values[m] = values[start];
      values[start] = tmp1;
    }
  }

  private double pvs(int depth, double alpha, double beta, int color) {
    double entry[] = {0, 0, -1, 0, 0}; // [value, depth, index, lower, upper]
    // only use a cached value if it was computed with the same depth - and
    // therefore the same ammount of knowledge (as this call) could be extracted.
    if (cache.containsKey(board.zobrist()) && (entry = cache.get(board.zobrist()))[1] >= depth) {
      cacheHits++;
      if (entry[3] >= beta) return entry[3];
      if (entry[4] <= alpha) return entry[4];
      alpha = Math.max(alpha, entry[3]);
      beta = Math.min(beta, entry[4]);
    } else cacheMisses++;
    MNKCell bc = null;
    double a = alpha, b = beta;

    minimaxed++;
    if (board.gameState() != MNKGameState.OPEN || depth == 0)
      // take into account the depth of the move when computing the score
      // expression:
      // color (1 for ME, -1 for ENEMY) * value (/ for ME, * for ENEMY) relativeDepth
      // a = color*evaluate() * (color < 0 ? 1d/Math.max(maxDepth-depth, 1) :
      // Math.max(maxDepth-depth, 1));
      a = evaluate(maxDepth - depth);
    else if (shouldHalt()) a = HALT;
    else if (board.getMarkedCells().length >= (K - 1) * 2
        && ((bc = findOneMoveWin(MY_WIN)) != null || (bc = findOneMoveLoss(ENEMY_WIN)) != null)) {
      board.markCell(bc.i, bc.j);
      a = -pvs(depth - 1, -beta, -alpha, -color);
      board.unmarkCell();
    } else {
      // {{{ sorting based on previous evaluations
      MNKCell[] cells = board.getFreeCells();
      double[] ratings = new double[cells.length];
      for (int i = 0; i < cells.length; i++) {
        long hash = board.nextZobrist(cells[i].i, cells[i].j);
        if (cache.containsKey(hash)) {
          double values[] = cache.get(hash);
          ratings[i] = values[1] >= depth ? -color * values[0] : 0;
          if (ratings[i] != 0) cacheHits++;
        } else {
          ratings[i] = 0;
          cacheMisses++;
        }
      }
      // }}}

      // {{{ minimax with alpha beta and pvs
      int i = 0;
      while (i < cells.length) {
        selectionSort(cells, ratings, i, color);
        MNKCell c = cells[i];
        board.markCell(c.i, c.j);
        double t;
        if (i > 0) {
          t = -pvs(depth - 1, -a - 1, -a, -color);
          if (t > a && t < b) t = -pvs(depth - 1, -b, -t, -color);
        } else t = -pvs(depth - 1, -b, -a, -color);
        board.unmarkCell();
        if (t == HALT) {
          if (a == alpha) return HALT;
          else break;
        }
        if (t > a) { // a = max(a,t)
          a = t;
          bc = c;
        }
        if (a >= b) {
          cutoff += board.getFreeCells().length - i;
          break;
        }
        b = a + 1; // null window for the next call
        i++;
      }
      // }}}
    }

    double lower = entry[3], upper = entry[4];
    if (a <= alpha) upper = a;
    if (a > alpha && a < beta) upper = lower = a;
    if (a >= beta) lower = a;

    // dummy values ignored regardless. Only the root value is taken into account
    if (bc == null && board.getFreeCells().length > 0) bc = board.getFreeCells()[0];

    double[] val = {a, depth, bc.i * minMN + bc.j, lower, upper};
    cache.put(board.zobrist(), val);
    return a;
  }

  public Pair<Double, int[]> iterativeDeepening() {
    int len = board.getFreeCells().length;
    double value = MIN;
    int cell[] = {0, 0};

    maxDepth = 1;
    while (!shouldHalt() && maxDepth < len) {
      double latest = pvs(maxDepth, MIN, MAX, 1);
      if (Math.abs(latest) == HALT) break;

      value = latest;
      int i = (int) cache.get(board.zobrist())[2];
      cell[0] = i / minMN;
      cell[1] = i % minMN;

      // stop the search if we found a certain win
      if (value >= winCutoff) break;

      // TODO: remove in production
      if (verbose)
        System.out.println(
            "minimax went to depth "
                + maxDepth
                + " with value: ("
                + value
                + ", ("
                + cell[0]
                + ","
                + cell[1]
                + "))");

      maxDepth++;
    }

    // TODO: remove in production
    if (maxDepth - 1 > board.getFreeCells().length) {
      throw new Error(
          "FATAL: iterativeDeepening exceeded allowable depth with a depth of: " + maxDepth);
    }

    return new Pair<>(value, cell);
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    startTime = System.currentTimeMillis();
    minimaxed = cutoff = seriesFound = evaluated = cacheHits = cacheMisses = 0;

    if (MC.length > 0)
      board.markCell(
          MC[MC.length - 1].i, MC[MC.length - 1].j); // keep track of the opponent's marks

    try {
      Pair<Double, int[]> result = iterativeDeepening();
      System.out.println(
          playerName()
              + "\t: visited "
              + minimaxed
              + " nodes, ended with result: ("
              + result.first
              + ", ("
              + result.second[0]
              + ","
              + result.second[1]
              + "))");
      System.out.println(
          playerName()
              + "\t: cut off "
              + cutoff
              + " branches ("
              + (cutoff / Math.max(minimaxed, 1) * 100)
              + "%)");
      System.out.println(
          playerName()
              + "\t: found a total of "
              + seriesFound
              + " free cells in series (up to k-1)");
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

      board.markCell(result.second[0], result.second[1]);
      return new MNKCell(result.second[0], result.second[1], MNKCellState.FREE);
    } catch (Exception e) {
      e.printStackTrace();
      return FC[new Random().nextInt(FC.length)];
    }
  }

  public String playerName() {
    return "Galileo Galilei";
  }
}
