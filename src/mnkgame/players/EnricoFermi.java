package mnkgame.players;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import mnkgame.*;

/*
 * Minimax + Alpha Beta + k-1,k-2,k-3 heuristics + euristic evaluation at fixed-depth
 * + caching (memoization)
 */
public class EnricoFermi implements MNKPlayer {
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
    public Board(int M, int N, int K) {
      super(M, N, K);
    }

    @Override
    public int hashCode() {
      // deep hash code of the contents of the array as it's multi dimensional,
      // otherwhise the hashes of the columns would not have been computed correctly
      return Arrays.deepHashCode(this.B);
    }

    @Override
    public String toString() {
      String str = "";
      for(int i = 0; i < M; i++) {
        for(int j = 0; j < N; j++)
          str += " " + (B[i][j] == MNKCellState.P1 ? 'x' : (B[i][j] == MNKCellState.P2 ? 'o' : '-'));
        str += '\n';
      }
      return str;
    }
  }
  // }}}

  // {{{ streak search
  private static enum StreakCellState {
    NONE(0), K3(1), K2(2), K1(3);
    private final int n;
    private StreakCellState(int n) {
      this.n = n;
    }

    public boolean isBetter(StreakCellState s) {
      return s.n > this.n;
    }
  }
  private static enum StreakDirection { HORIZONTAL, VERTICAL, DIAGONAL, ANTIDIAGONAL; }

  private static int streak = 0, streak_free = 0;
  private static MNKCellState streak_enemy;
  private static StreakCellState[] streak_states = StreakCellState.values();
  private static StreakCellState[][] SB;
  private static StreakDirection streakDirection;

  // adds the cell at the given position to the curret streak and computes if
  // we should care about it. Returns true if the streak has been interrupted
  private boolean checkStreakCell(int i, int j) {
    MNKCellState c = board.cellState(i, j);
    if(c == streak_enemy)
      return true;
    else {
      streak++;
      if(c == MNKCellState.FREE)
        streak_free++;
    }

    if(streak == K && streak_free <= 3)
      markStreakCell(i, j);
    return false;
  }

  private void markStreakCell(int i, int j) {
    StreakCellState s = streak_states[streak_free];
    switch(streakDirection) {
      case HORIZONTAL:
        for(int x = i; x > i-K; x--)
          if(SB[x][j].isBetter(s))
            SB[x][j] = s;
      break;

      default:
      break;
    }

    // decrease the streak by 1 so we count further cells;
    streak--;
    streak_free--;
  }

  private void resetStreak() {
    streak = streak_free = 0;
  }

  private void findStreaks(MNKCellState player) {
    streak_enemy = opponent(player);

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        SB[i][j] = StreakCellState.NONE;

    // iterate over all rows
    for(int i = 0; i < M; i++) {
      resetStreak();
      for(int j = 0; j < N; j++)
        if(checkStreakCell(i, j))
          resetStreak();
    }
  }
  // }}}

  // {{{ series
  private class Series {
    private static class CellState {
      public int series; // 1-6 int to mark the current best series for the given cell
      public MNKCell cell;

      public CellState(int series, MNKCell cell) {
        this.series = series;
        this.cell = cell;
      }

      @Override
      public String toString() {
        return "CellState[series=" + series + ",cell=" + cell + "]";
      }
    }

    public static class Result {
      public List<MNKCell> k1;
      public List<MNKCell> kn;

      public Result(List<MNKCell> k1, List<MNKCell> kn) {
        this.k1 = k1;
        this.kn = kn;
      }
    }

    private static Board board = null;
    private static int M, N, K;
    private static MNKCellState player = null, opponent = null;
    private static CellState[][] cells;

    private static List<List<MNKCell>> series =
        new ArrayList<List<MNKCell>>(); // k-(1..3) streak cells for the current player/opponent
    // k-n where n is odd means k-n streak for the current player
    // k-n where n is even means k-(n/2) streak for the opponent

    private static MNKCellState curr_player = MNKCellState.FREE;
    private static List<MNKCell> curr_streak = new ArrayList<>();
    private static List<MNKCell> curr_streak_free = new ArrayList<>();

    private static void reset() {
      curr_player = MNKCellState.FREE;
      curr_streak.clear();
      curr_streak_free.clear();
    }

    private static void addCell(int i, int j, int s) {
      // first check if we have to remove the cell
      CellState prev;
      if ((prev = cells[i][j]) != null && prev.series > s)
        // TODO: improve cost: MxN matrix with double-linked-list structure
        series.get(prev.series).remove(prev.cell);

      MNKCell c = new MNKCell(i, j, MNKCellState.FREE);
      cells[i][j] = new CellState(s, c);
      series.get(s).add(c);
    }

    private static int index(int streak, MNKCellState player) {
      if (streak == 1) return player == opponent ? 1 : 0;
      else if (streak == 2) return 2 + (player == opponent ? 1 : 0);
      else return 4 + (player == opponent ? 1 : 0);
    }

    private static void checkCell(int i, int j) {
      MNKCellState s = board.cellState(i, j);
      MNKCell c = new MNKCell(i, j, s);
      // when we have no current player we wanna fill that asap
      // take this board as an example(k=3):
      // - - - -  | the first two cells in the first row will be checked and
      // - - - -  | added, later we find the x and want to make sure the previously
      // x - - -  | discorvered cells are used in the x player streak.
      // o - - -  |
      if (curr_player == MNKCellState.FREE && s != MNKCellState.FREE) curr_player = s;

      if (s == curr_player || s == MNKCellState.FREE) {
        // possible streak cell
        curr_streak.add(c);
        if (s == MNKCellState.FREE) curr_streak_free.add(c);
      } else {
        // opponent marked cell swap streak data
        curr_streak.clear();
        curr_streak_free.clear();
        curr_streak.add(c);
        curr_player = s;
      }

      // we have a k-(n<=3) streak assigned to a player
      if (curr_streak.size() >= K
          && curr_streak_free.size() <= 3
          && curr_player != MNKCellState.FREE) {
        for (MNKCell cc : curr_streak_free)
          addCell(cc.i, cc.j, index(curr_streak_free.size(), curr_player));
      }
      // remove the last free cell as it'll be useless in the next iteration,
      // as we want a k-3 streak at min
      if (curr_streak_free.size() >= 3) {
        // ignore the first free cell in the streak to allow the lookup of other
        // possible streaks. Take this board as an example(k=4):
        // - - - - o - | in the first 4 iterations we find 4 free cells but the
        // - - - - o - | curr_player is not clear so we remove the (1,1) cell and
        // - - - - - - | move onto the 5th iteration, where we find a marked cell
        // - - - - - - | and have a streak match.
        // x - - - - - |
        // x - - - - - |
        // Another example: (which shows why we need that loop)
        // TODO: cleanup comments
        // x           (k=5)
        // -
        // -
        // -
        // -
        // x
        // x
        MNKCell r = curr_streak_free.remove(0);
        while (curr_streak.get(0) != r) curr_streak.remove(0);
        curr_streak.remove(0);
      }
    }

    // returns a pair of (k-1, k-n) cells. We separate the k-1 cells (for both players)
    // as these are the ones we can mark instantly without developing the whole
    // minimax tree as we know they are the possible available moves (either give
    // a win or prevent a loss)
    public static Result findSeries(Board b, MNKCellState p) {
      board = b;
      M = board.M;
      N = board.N;
      K = board.K;
      cells = new CellState[M][N];
      series.clear();
      for (int i = 0; i < 6; i++) series.add(new ArrayList<MNKCell>());

      player = p;
      opponent = opponent(p);

      // iterate over all rows
      for (int i = 0; i < M; i++) {
        reset();
        for (int j = 0; j < N; j++) checkCell(i, j);
      }

      // iterate over all columns
      for (int j = 0; j < N; j++) {
        reset();
        for (int i = 0; i < M; i++) checkCell(i, j);
      }

      // iterate over all diagonals
      int nDiagonals = (Math.min(N, M) - K) * 2 + 1;
      for (int x = 0; x < nDiagonals; x++) {
        reset();
        int i = 0, j = 0;
        if (x != 0 && x % 2 == 0) i = x / 2;
        else if (x != 0) j = x;

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

      // k[0] and k[1] are k-1 series for both player, same relevance but
      // we wanna keep the k-1 for urselves first as they lead to win and
      // not loss defense
      List<MNKCell> k1 = new ArrayList<>(series.get(0));
      k1.addAll(series.get(1));

      // parepare the results
      List<MNKCell> kn = new LinkedList<>();
      for (int i = 2; i < 6; i++) kn.addAll(series.get(i));

      // increase the amount of series found here as later we'll also add
      // non-series cells to kn
      series_found += k1.size() + kn.size();

      for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
          if (cells[i][j] == null && b.cellState(i, j) == MNKCellState.FREE)
            kn.add(new MNKCell(i, j, MNKCellState.FREE));

      return new Result(k1, kn);
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

    private static int curr_streak, curr_streak_free;

    private static double checkCell(int i, int j) {
      MNKCellState c = board.cellState(i, j);
      if (c == player || c == MNKCellState.FREE) {
        curr_streak++;
        curr_streak_free++;
      } else curr_streak = curr_streak_free = 0;

      double value = 0;
      if (curr_streak >= K) {
        // we know for a fact that curr_stack_free > 0, oterwhise this method
        // won't be called as the evaluation can be done in a deterministic way
        value = Math.log(curr_streak / curr_streak_free) / klog;
        // should be 1 when curr_streak_free == 1 which is the best value and
        // the greatest chance to win (only need to mark 1 cell) whereas
        // should get closer to 0 as curr_streak_free increases. It's effectively
        // log_k(curr_streak/curr_streak_free)
        curr_streak--;
        curr_streak_free--;
      }
      return value;
    }

    private static void reset() {
      curr_streak = curr_streak_free = 0;
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
      for (int i = 0; i < M && M - i + curr_streak >= K; i++) {
        reset();
        for (int j = 0; j < N && N - j + curr_streak >= K; j++) chance += checkCell(i, j);
      }

      // check all rows
      for (int j = 0; j < N && N - j + curr_streak >= K; j++) {
        reset();
        for (int i = 0; i < M && M - i + curr_streak >= K; i++) chance += checkCell(i, j);
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

  private final int HALT = 1337;
  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, ENEMY_WIN;
  private int M, N, K;
  private long start_time, timeout;
  private int maxDepth;
  private Board board;
  private HashMap<String, Pair<Integer, Integer>> cache = new HashMap<>();

  // NOTE: profiling
  private static int visited, series_found, evaluated, cache_hits, cache_misses;
  private static boolean clear = false;

  enum Action {
    MINIMIZE,
    MAXIMIZE
  }

  private static Action opposite(final Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  private static MNKCellState opponent(final MNKCellState p) {
    return p == MNKCellState.P1 ? MNKCellState.P2 : MNKCellState.P1;
  }

  private boolean shouldHalt() {
    // TODO: tweak values
    return (System.currentTimeMillis() - start_time) / 1000.0 > timeout * 0.95;
  }

  // finds the first cell needed to copmlete a K-1 streak in any possible direction
  private MNKCell findOneMoveWin(final MNKGameState win_state) {
    for (MNKCell c : board.getFreeCells()) {
      MNKGameState result = board.markCell(c.i, c.j);
      board.unmarkCell();
      if (result == win_state) return c;
    }
    return null;
  }

  private MNKCell pickRandomNonClosingCell(MNKCell previous) {
    for (MNKCell c : board.getFreeCells()) {
      // avoiding the should_halt check here, kinda superflous
      MNKGameState result = board.markCell(c.i, c.j);
      board.unmarkCell();
      if (result == MNKGameState.OPEN
          && (previous == null || previous.i != c.i || previous.j != c.j)) return c;
    }
    return null;
  }

  // finds the first cell the enemy needs to copmlete a K-1 streak in any possible direction
  private MNKCell findOneMoveLoss(final MNKGameState loss_state) {
    MNKCell random_cell = null;
    if (board.getFreeCells().length == 1
        || (random_cell = pickRandomNonClosingCell(null)) == null)
      return null; // cannot check for enemy's next move when it doesn't exist

    board.markCell(random_cell.i, random_cell.j);
    MNKCell c = findOneMoveWin(loss_state);
    board.unmarkCell(); // remove the marked random_cell
    if (c != null) return c;

    // test the random_cell we selected at first. It may be a one-move loss cell
    // get a new random cell different from the previous and call it safe_cell
    MNKCell cc = pickRandomNonClosingCell(random_cell);
    if (board.markCell(cc.i, cc.j) != MNKGameState.OPEN) {
      // random_cell puts us in a draw, ignore that
      board.unmarkCell();
      return null;
    }
    MNKGameState result = board.markCell(random_cell.i, random_cell.j); // let the enemy take the random ane
    board.unmarkCell();
    board.unmarkCell();
    return result == loss_state ? random_cell : null;
  }

  private static final int winCutoff = 2; // 2 represents a certain win, and therefore we can cutoff search
  // evaluate a state (either heuristically or in a deterministic way) regardless
  // of its depth. The depth will be taken into account later to allow for caching
  private int evaluate(final Board board) {
    if (board.gameState() == MY_WIN) return 2;
    else if (board.gameState() == ENEMY_WIN) return -2;
    else if (board.gameState() == MNKGameState.DRAW) return 0;
    else {
      evaluated++;
      // keep the heuristic evaluation between 1 and -1
      return (int) Math.min(Math.max(Chances.winningChances(board, ME) - Chances.winningChances(board, ENEMY), 1), -1);
    }
  }

  private Tuple<Integer, Integer, MNKCell> mem(Tuple<Integer, Integer, MNKCell> result) {
    if(result.second > 1)
      cache.put(board.toString(), new Pair<>(result.first.intValue(), board.getMarkedCells().length));
    return result;
  }

  private Pair<Integer, Integer> unmem() {
    if(cache.containsKey(board.toString()))
      return cache.get(board.toString());
    else
      return null;
  }

  private Tuple<Integer, Integer, MNKCell> minimax(Action action, int depth, double a, double b) {
    Pair<Integer, Integer> cached;
    if (depth > 0 && (cached = unmem()) != null) {
      cache_hits++;
      return new Tuple<>(
          cached.first, cached.second, null);
    } else cache_misses++;

    if (board.gameState() != MNKGameState.OPEN || depth == maxDepth)
      // return the evaluation of the current board with the last marked cell
      // as the decision which has brough up to this game state
      return mem(
          new Tuple<>(
              evaluate(board), depth, board.getMarkedCells()[board.getMarkedCells().length - 1]));

    if (shouldHalt()) return new Tuple<>(HALT, depth, board.getFreeCells()[0]);

    visited++;

    MNKCell omc = null; // one move win/loss cell
    if (board.getMarkedCells().length >= (K - 1) * 2
        && // only check for one-win-moves
        // if there have been placed
        // enough cells to make one happen
        ((omc = findOneMoveWin(action == Action.MAXIMIZE ? MY_WIN : ENEMY_WIN)) != null
            || (omc = findOneMoveLoss(action == Action.MAXIMIZE ? ENEMY_WIN : MY_WIN))
                != null)) {
      // if we know we are acting on the root we can just return without
      // evaluating the value of the move, as this won't be used anywhere
      series_found++; // found once cell in a k-1 series

      board.markCell(omc.i, omc.j);
      Tuple<Integer, Integer, MNKCell> result = minimax(opposite(action), depth + 1, a, b);
      board.unmarkCell();
      return mem(new Tuple<>(result.first, result.second, omc));
    }

    double best = action == Action.MAXIMIZE ? -Double.MAX_VALUE : Double.MAX_VALUE;
    int best_value = Integer.MIN_VALUE, best_depth = depth + 1;
    MNKCell best_cell = null;
    // TODO: (?) shuffle(series.noseries);
    for (MNKCell c : board.getFreeCells()) {
      board.markCell(c.i, c.j);
      Tuple<Integer, Integer, MNKCell> result = minimax(opposite(action), depth + 1, a, b);
      board.unmarkCell();

      if(result.first == HALT) break;

      double value =
          action == Action.MAXIMIZE ? (double) result.first / result.second : (double) result.first * result.second;

      if ((action == Action.MAXIMIZE && value > best)
          || (action == Action.MINIMIZE && value < best)) {
        best = value;
        best_value = result.first;
        best_depth = result.second;
        best_cell = c;
      }
      if (action == Action.MAXIMIZE && best > a) a = best;
      else if (action == Action.MINIMIZE && best < b) b = best;

      if (b <= a) break;
    }
    return mem(new Tuple<>(best_value, best_depth, best_cell));
  }

  public void initPlayer(int M, int N, int K, boolean first, int timeout_in_secs) {
    this.M = M;
    this.N = N;
    this.K = K;
    this.timeout = timeout_in_secs;

    MY_WIN = first ? MNKGameState.WINP1 : MNKGameState.WINP2;
    ENEMY_WIN = first ? MNKGameState.WINP2 : MNKGameState.WINP1;
    ME = first ? MNKCellState.P1 : MNKCellState.P2;
    ENEMY = first ? MNKCellState.P2 : MNKCellState.P1;
    board = new Board(M, N, K);
    cache.clear();

    // clear screen
    if(clear) {
      System.out.print("\033[H\033[2J");  
      System.out.flush();
    }
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    visited = series_found = evaluated = cache_hits = cache_misses = 0;
    maxDepth = 6; // TODO: dynamic

    if (MC.length > 0)
      board.markCell(
          MC[MC.length - 1].i, MC[MC.length - 1].j); // keep track of the opponent's marks

    try {
      Tuple<Integer, Integer, MNKCell> result =
          minimax(Action.MAXIMIZE, 0, -Double.MAX_VALUE, Double.MAX_VALUE);
      System.out.println(
          playerName() + "\t: visited " + visited + " nodes, ended with result: " + result);
      System.out.println(
          playerName()
              + "\t: found a total of "
              + series_found
              + " free cells in series (up to k-3)");
      System.out.println(playerName() + "\t: heuristically evaluated " + evaluated + " boards");
      int perc = (int) (((double) cache_hits / (cache_misses + cache_hits)) * 100);
      System.out.println(
          playerName()
              + "\t: cached "
              + cache.size()
              + " elements, hit: "
              + cache_hits
              + ", misses: "
              + cache_misses
              + ". rate: "
              + perc
              + "%");
      System.out.println(
          playerName()
              + "\t: took " + (System.currentTimeMillis()-start_time)/1000d + "s"
          );

      // TODO: remove in prod
      if (FC.length != board.getFreeCells().length) {
        System.out.println("FATAL: minimax didn't clean the board");
        return FC[0];
      }

      board.markCell(result.third.i, result.third.j);
      return result.third;
    } catch (Exception e) {
      e.printStackTrace();
      return FC[new Random().nextInt(FC.length)];
    }
  }

  public String playerName() {
    return "Enrico Fermi";
  }
}
