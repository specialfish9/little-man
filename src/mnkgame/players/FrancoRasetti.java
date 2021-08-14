package mnkgame.players;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import mnkgame.*;

/*
 * Minimax + Alpha Beta + Iterative Deepening + k-1 heuristics
 * + euristic evaluation at iterative depth + caching (memoization)
 */
public class FrancoRasetti implements MNKPlayer {
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
  private boolean finite; // when we have reached max depth
  private Board board;

  private HashMap<String, Pair<Integer, Integer>> cache = new HashMap<>();

  // NOTE: profiling
  private static int visited, series_found, evaluated, cache_hits, cache_misses;
  private static boolean clear = false;

  public void initPlayer(int M, int N, int K, boolean first, int timeoutInSecs) {
    this.M = M;
    this.N = N;
    this.K = K;
    this.timeout = timeoutInSecs;

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

  enum Action {
    MINIMIZE,
    MAXIMIZE
  }

  private static Action opposite(final Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  private boolean shouldHalt() {
    // TODO: tweak values
    return (System.currentTimeMillis() - start_time) / 1000.0 > timeout * 0.98; // livin' on the edge
  }

  private boolean shouldHalt(long endTime) {
    return System.currentTimeMillis()>=endTime;
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

  // evaluate a state (either heuristically or in a deterministic way) regardless
  // of its depth. The depth will be taken into account later to allow for caching
  private static final int winCutoff = 2; // 2 represents a certain win, and therefore we can cutoff search
  // values for board evaluation:
  // - negative for heuristics
  // - positive for deterministic
  // - (deterministic) 1 = loss, 2 = draw, 3 = win
  // - (heuristic) -1 = loss, -2 = draw, -3 = win
  // we use the oddness of the number and its sign to slightly manipulate it and
  // use it for alpha/beta values
  private int evaluate(final Board board) {
    // check for draws first, most lickely
    if (board.gameState() == MNKGameState.DRAW) return 2;
    else if (board.gameState() == MY_WIN) return 3;
    else if (board.gameState() == ENEMY_WIN) return 1;
    else {
      // TODO: only check caching here. FIX?
      if (cache.containsKey(board.toString())) {
        cache_hits++;
        return cache.get(board.toString()).first;
      } else cache_misses++;

      evaluated++;
      // keep the heuristic evaluation between 1 and -1
      int res = (int) Math.min(Math.max(Chances.winningChances(board, ME) - Chances.winningChances(board, ENEMY), -1), 1);
      // System.out.println("user points: " + Chances.winningChances(board, ME));
      // System.out.println("enemy points: " + Chances.winningChances(board, ENEMY));
      cache.put(board.toString(), new Pair<>(res, 13377)); // ultra fucked up
      return res;
    }
  }

  private Pair<Integer, Integer> minimax(Action action, int depth, long endTime, double a, double b) {
    if (board.gameState() != MNKGameState.OPEN || depth == 0) {
      if(depth != 0)
        finite = true;
      // return the evaluation of the current board with the last marked cell
      // as the decision which has brough up to this game state
      return new Pair<>(evaluate(board), depth);
    }

    if (shouldHalt(endTime)) return new Pair<>(HALT, depth);

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
      Pair<Integer, Integer> result = minimax(opposite(action), depth - 1, endTime, a, b);
      board.unmarkCell();
      return new Pair<>(result.first, result.second);
    }

    double best = action == Action.MAXIMIZE ? -Double.MAX_VALUE : Double.MAX_VALUE;
    int best_value = Integer.MIN_VALUE, best_depth = depth - 1;
    int len = board.getFreeCells().length;
    // TODO: (?) shuffle(series.noseries);
    for (MNKCell c : board.getFreeCells()) {
      board.markCell(c.i, c.j);
      Pair<Integer, Integer> result = minimax(opposite(action), depth - 1, endTime, a, b);
      board.unmarkCell();

      if(result.first == HALT) break;

      // refer to the comment on the evaluate function
      // TODO: factor in len-depth
  // values for board evaluation:
  // - negative for heuristics
  // - positive for deterministic
  // - (deterministic) 1 = loss, 2 = draw, 3 = win
  // - (heuristic) -1 = loss, -2 = draw, -3 = win
  // we use the oddness of the number and its sign to slightly manipulate it and
  // use it for alpha/beta values
      int value = result.first % 2 == 0
        ? 0 /* draw */
        : (Math.abs(result.first) == 1 ? - );
          // action == Action.MAXIMIZE ? (double) result.first / (len-result.second) : (double) result.first * (len-result.second);
      // used to distinguish actual draw from heuristic draw
      boolean heuristic = result.first % 2 == 0;
      if(!heuristic)
        value *= 2; // certain win/loss = 2 or -2, heuristic = 1 or -1


      if ((action == Action.MAXIMIZE && value > best)
          || (action == Action.MINIMIZE && value < best)) {
        best = value;
        best_value = result.first;
        best_depth = result.second;
      }
      if (action == Action.MAXIMIZE && best > a) a = best;
      else if (action == Action.MINIMIZE && best < b) b = best;

      // TODO: check back
      if (b <= a /*|| (action == Action.MAXIMIZE && best_value >= winCutoff)*/) break;
    }
    return new Pair<>(best_value, best_depth);
  }

  public int iterativeDeepeningRating(long endTime) {
    int depth = 1;
    finite = false;
    int result = Integer.MIN_VALUE;

    // as long as we can deepen the tree for this move
    while(!shouldHalt(endTime) && !finite) {
      // TODO: remove in production
      if(depth-1 > board.getFreeCells().length)
        System.out.println("FATAL: iterativeDeepening exceeded allowable depth with a depth of: " + depth);

      Pair<Integer, Integer> val = minimax(Action.MINIMIZE, depth, endTime, -Double.MAX_VALUE, Double.MAX_VALUE);
      System.out.println("minimax at depth " + depth + " has value: " + val);

      if(val.first == HALT)
        break;
      else if(val.first >= winCutoff)
        return val.first;
      else
        result = val.first; // always keep the last result as it's the most accurate

      depth++;
    }
    return result;
  }

  public MNKCell find() {
    int len = board.getFreeCells().length;
    int best_value = Integer.MIN_VALUE;
    MNKCell best_cell = null;
    for(MNKCell c : board.getFreeCells()) {
      long endTime = System.currentTimeMillis()+(timeout*980)/len;
      board.markCell(c.i, c.j);
      int score = iterativeDeepeningRating(endTime);

      board.unmarkCell();

      if(score >= winCutoff)
        return c;

      if(score > best_value) {
        best_value = score;
        best_cell = c;
      }
    }
    return best_cell;
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    visited = series_found = evaluated = cache_hits = cache_misses = 0;

    if (MC.length > 0)
      board.markCell(
          MC[MC.length - 1].i, MC[MC.length - 1].j); // keep track of the opponent's marks

    try {
      MNKCell result = find();
      System.out.println(
          playerName() + "\t: visited " + visited + " nodes, ended with result: (x,x," + result + ")");
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
