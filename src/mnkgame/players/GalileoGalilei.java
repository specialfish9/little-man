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

    public Board(int M, int N, int K, int minMN, MNKCellState me) {
      super(M, N, K);
      this.minMN = minMN;
      this.me = me;
    }

    @Override
    public MNKGameState markCell(int i, int j) {
      // mind the order of the calls
      key = nextZobrist(i, j);
      previousValues.push(value);
      value -= eval(i, j);
      MNKGameState result = super.markCell(i, j);
      value += eval(i, j);
      return result;
    }

    @Override
    public void unmarkCell() {
      // mind the order of the calls
      MNKCell last = MC.getLast();
      super.unmarkCell();
      key = nextZobrist(last.i, last.j);
      value = previousValues.pop();
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

    private double pushCell(final MNKCellState state) {
      if (queue.size() >= K) // useless >
      popCell();

      if (state == MNKCellState.FREE) queueFree++;
      else if (state == MNKCellState.P1) queueP1++;
      else if (state == MNKCellState.P2) queueP2++;
      queue.add(state);
      int sign = me == MNKCellState.P1 ? 1 : -1;
      if (queueP1 + queueFree == K) return sign * (seriesValue(queueFree) + (queueP1 * 1d / K));
      else if (queueP2 + queueFree == K)
        return -sign * (seriesValue(queueFree) + (queueP2 * 1d / K));
      else return 0;
    }

    private int seriesValue(final int free) {
      if (K > 4 && free == 3) return 1000;
      if (K > 3 && free == 2) return 10000;
      if (K > 2 && free == 1) return 100000;

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

  private static final double HALT = Double.MIN_VALUE;
  private static final double MIN = -Double.MAX_VALUE, MAX = Double.MAX_VALUE;
  private static double winCutoff; // represents a certain win, and therefore we can cutoff search

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
  private MNKCell bestCell = null;

  // NOTE: profiling
  private static int minimaxed,
      failed,
      failedDepth,
      cutoff,
      seriesFound,
      evaluated,
      cacheHits,
      cacheMisses;
  private static boolean clear = false, verbose = true;

  public void initPlayer(int M, int N, int K, boolean first, int timeoutInSecs) {
    this.M = M;
    this.N = N;
    this.K = K;
    this.minMN = Math.min(M, N);
    this.winCutoff = MAX / (M * N);
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

    board.markCell(6, 4);
    System.out.println(board + " -> " + board.value());
    board.markCell(6, 5);
    System.out.println(board + " -> " + board.value());
    board.markCell(5, 4);
    System.out.println(board + " -> " + board.value());
    board.markCell(4, 4);
    System.out.println(board + " -> " + board.value());

    board.markCell(6, 2);
    System.out.println(board + " -> " + board.value());
    board.unmarkCell();

    board.markCell(5, 3);
    System.out.println(board + " -> " + board.value());
    board.unmarkCell();

    board.unmarkCell();
    board.unmarkCell();
    board.unmarkCell();
    board.unmarkCell();
  }

  private boolean shouldHalt() {
    // TODO: tweak values
    return (System.currentTimeMillis() - startTime) / 1000.0 >= timeout * 0.9; // livin' on the edge
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
    board.markCell(cc.i, cc.j);
    MNKGameState result =
        board.markCell(randomCell.i, randomCell.j); // let the enemy take the random ane
    board.unmarkCell();
    board.unmarkCell();
    return result == lossState ? randomCell : null;
  }

  // evaluate a state (either heuristically or in a deterministic way) regardless
  // of its depth. The depth is also taken into account
  private double evaluate() {
    // check for draws first, most lickely
    double gameDepth = Math.max(board.getMarkedCells().length, 1);
    MNKGameState state = board.gameState();
    if (state == MNKGameState.DRAW) return 0;
    else if (state == MY_WIN) return (MAX - 1) / gameDepth;
    else if (state == ENEMY_WIN) return -(MIN + 1) / gameDepth;
    else {
      evaluated++;
      return board.value() / gameDepth;
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

  private double pvs(int searchDepth, double alpha, double beta, int color) {
    double entry[] = {0, 0, -1, 0, 0}; // [value, depth, index, lower, upper]
    // only use a cached value if it was computed with the same depth - and
    // therefore the same ammount of knowledge (as this call) could be extracted.
    if (cache.containsKey(board.zobrist())
        && (entry = cache.get(board.zobrist()))[1] >= searchDepth) {
      cacheHits++;
      if (entry[3] >= beta) return entry[3];
      if (entry[4] <= alpha) return entry[4];
      alpha = Math.max(alpha, entry[3]);
      beta = Math.min(beta, entry[4]);
    } else cacheMisses++;
    MNKCell bc = null;
    double a = alpha, b = beta;

    minimaxed++;
    if (board.gameState() != MNKGameState.OPEN || searchDepth <= 0) {
      a = color * evaluate();
    } else if (shouldHalt()) a = HALT;
    else if (board.getMarkedCells().length >= (K - 1) * 2 - 1
        && ((bc = findOneMoveWin(MY_WIN)) != null || (bc = findOneMoveLoss(ENEMY_WIN)) != null)) {
      seriesFound++;
      board.markCell(bc.i, bc.j);
      a = -pvs(searchDepth - 1, -beta, -alpha, -color);
      board.unmarkCell();
    } else {
      // {{{ sorting based on previous evaluations
      MNKCell[] cells = board.getFreeCells();
      double[] ratings = new double[cells.length];
      for (int i = 0; i < cells.length; i++) {
        long hash = board.nextZobrist(cells[i].i, cells[i].j);
        if (cache.containsKey(hash)) {
          double values[] = cache.get(hash);
          ratings[i] = values[1] >= searchDepth ? values[0] : 0;
        } else ratings[i] = 0;

        // TODO: remove in production
        if (ratings[i] != 0) cacheHits++;
        else cacheMisses++;
      }
      // }}}

      // {{{ minimax with alpha beta and pvs
      int i = 0;
      while (i < cells.length) {
        selectionSort(cells, ratings, i, color);
        MNKCell c = cells[i];
        board.markCell(c.i, c.j);
        double t;
        if (i == 0) t = -pvs(searchDepth - 1, -b, -a, -color);
        else {
          t = -pvs(searchDepth - 1, -a - 1, -a, -color); // null-window search
          if (t > a && t < b) {
            failed++;
            failedDepth += searchDepth - 1;
            t =
                -pvs(
                    searchDepth - 1,
                    -b,
                    -t,
                    -color); // if we fail between a and b do a proper search
          }
        }
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
        i++;
      }
      // }}}
    }
    if (a == HALT) return HALT;

    double lower = entry[3], upper = entry[4];
    if (a <= alpha) upper = a;
    if (a > alpha && a < beta) upper = lower = a;
    if (a >= beta) lower = a;

    double[] val = {a, searchDepth, bc == null ? -1 : bc.i * minMN + bc.j, lower, upper};
    cache.put(board.zobrist(), val);
    bestCell = bc;
    return a;
  }

  public Pair<Double, MNKCell> iterativeDeepening() {
    int len = board.getFreeCells().length;
    double value = MIN;
    MNKCell bc = null;

    maxDepth = 1;
    while (!shouldHalt() && maxDepth <= len) {
      double latest = pvs(maxDepth, MIN, MAX, 1);
      if (Math.abs(latest) == HALT) break;

      value = latest;
      bc = bestCell;

      // stop the search if we found a certain win
      if (value >= winCutoff) break;

      // TODO: remove in production
      if (verbose)
        System.out.println(
            "minimax went to depth " + maxDepth + " with value: (" + value + ", " + bc + ")");

      maxDepth++;
    }

    // TODO: remove in production
    if (maxDepth - 1 > board.getFreeCells().length) {
      throw new Error(
          "FATAL: iterativeDeepening exceeded allowable depth with a depth of: " + maxDepth);
    }

    return new Pair<>(value, bc);
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    startTime = System.currentTimeMillis();
    minimaxed =
        failed = failedDepth = cutoff = seriesFound = evaluated = cacheHits = cacheMisses = 0;

    if (MC.length > 0)
      board.markCell(
          MC[MC.length - 1].i, MC[MC.length - 1].j); // keep track of the opponent's marks

    try {
      Pair<Double, MNKCell> result = iterativeDeepening();
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
          playerName() + "\t: took " + (System.currentTimeMillis() - startTime) / 1000d + "s");

      // TODO: remove in prod
      if (FC.length != board.getFreeCells().length) {
        System.out.println("FATAL: minimax didn't clean the board");
        return FC[0];
      }

      board.markCell(result.second.i, result.second.j);
      return result.second;
    } catch (Exception e) {
      e.printStackTrace();
      MNKCell c = FC[new Random().nextInt(FC.length)];
      board.markCell(c.i, c.j);
      return c;
    }
  }

  public String playerName() {
    return "Galileo Galilei";
  }
}
