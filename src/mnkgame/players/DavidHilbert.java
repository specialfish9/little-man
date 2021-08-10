package mnkgame.players;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import mnkgame.*;

/*
 * Minimax + Alpha Beta + one-win-euristics (next-move wind/loss interception)
 * + radious-checking cells + 
 */
public class DavidHilbert implements MNKPlayer {
  private static class Pair<A, B> {
    public A first;
    public B second;

    Pair(A first, B second) {
      this.first = first;
      this.second = second;
    }

    @Override
    public String toString() {
      return "(" + first + "," + second + ")";
    }
  }

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

    private static MNKBoard board = null;
    private static int M, N, K;
    private static MNKCellState player = null, opponent = null;
    private static CellState[][] cells;

    private static List<List<MNKCell>> series = new ArrayList<List<MNKCell>>(); // k-(1..3) streak cells for the current player/opponent 
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
      if((prev = cells[i][j]) != null && prev.series > s)
        // TODO: improve cost: MxN matrix with double-linked-list structure
        series.get(prev.series).remove(prev.cell);

      MNKCell c = new MNKCell(i, j, MNKCellState.FREE);
      cells[i][j] = new CellState(s, c);
      series.get(s).add(c);
    }

    private static int index(int streak, MNKCellState player) {
      if(streak == 1)
        return player == opponent ? 1 : 0;
      else if(streak == 2)
        return 2 + (player == opponent ? 1 : 0);
      else
        return 4 + (player == opponent ? 1 : 0); 
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
      if(curr_player == MNKCellState.FREE && s != MNKCellState.FREE)
        curr_player = s;

      if(s == curr_player || s == MNKCellState.FREE) {
        // possible streak cell
        curr_streak.add(c);
        if(s == MNKCellState.FREE)
          curr_streak_free.add(c);
      } else {
        // opponent marked cell swap streak data
        curr_streak.clear();
        curr_streak_free.clear();
        curr_streak.add(c);
        curr_player = s;
      }

      // we have a k-(n<=3) streak assigned to a player
      if(curr_streak.size() >= K && curr_streak_free.size() <= 3 &&
          curr_player != MNKCellState.FREE) {
        for(MNKCell cc : curr_streak_free)
          addCell(cc.i, cc.j, index(curr_streak_free.size(), curr_player));
      }
      // remove the last free cell as it'll be useless in the next iteration, 
      // as we want a k-3 streak at min
      if(curr_streak_free.size() >= 3) {
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
        while(curr_streak.get(0) != r)
          curr_streak.remove(0);
        curr_streak.remove(0);
      }
    }

    // returns a pair of (k-1, k-n) cells. We separate the k-1 cells (for both players)
    // as these are the ones we can mark instantly without developing the whole
    // minimax tree as we know they are the possible available moves (either give
    // a win or prevent a loss)
    public static Result findSeries(MNKBoard b, MNKCellState p) {
      board = b; 
      M = board.M; N = board.N; K = board.K;
      cells = new CellState[M][N];
      series.clear();
      for(int i = 0; i < 6; i++)
        series.add(new ArrayList<MNKCell>());

      player = p;
      opponent = opponent(p);

      // iterate over all rows
      for(int i = 0; i < M; i++) {
        reset();
        for(int j = 0; j < N; j++)
          checkCell(i, j);
      }
 
      // iterate over all columns
      for(int j = 0; j < N; j++) {
        reset();
        for(int i = 0; i < M; i++)
          checkCell(i, j);
      }

      // iterate over all diagonals
      int nDiagonals = (Math.min(N, M) - K)*2 + 1;
      for(int x = 0; x < nDiagonals; x++) {
        reset();
        int i = 0, j = 0;
        if(x != 0 && x % 2 == 0)
          i = x/2;
        else if (x != 0)
          j = x;

        while(i < M && j < N) {
          checkCell(i, j);
          j++; i++;
        }
      }

      // iterate over all counter diagonals
      for(int x = 0; x < nDiagonals; x++) {
        reset();
        int i = 0, j = N-1;
        if(x != 0 && x % 2 == 0)
          i = x/2;
        else if (x != 0)
          j = N-1-x;

        while(i < M && j >= 0) {
          checkCell(i, j);
          j--; i++;
        }
      }

      // k[0] and k[1] are k-1 series for both player, same relevance but
      // we wanna keep the k-1 for urselves first as they lead to win and
      // not loss defense
      List<MNKCell> k1 = new ArrayList<>(series.get(0));
      k1.addAll(series.get(1));

      // parepare the results
      List<MNKCell> kn = new LinkedList<>();
      for(int i = 2; i < 6; i++)
        kn.addAll(series.get(i));

      // increase the amount of series found here as later we'll also add 
      // non-series cells to kn
      series_found += k1.size() + kn.size();
    

      for(int i = 0; i < M; i++)
        for(int j = 0; j < N; j++)
          if(cells[i][j] == null && b.cellState(i, j) == MNKCellState.FREE)
            kn.add(new MNKCell(i, j, MNKCellState.FREE));

      return new Result(k1, kn);
    }
  }

  private final double RANK_CONSTANT = 10;
  private final double HALT = Double.MIN_VALUE;
  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, ENEMY_WIN;

  private int M, N, K;
  private long start_time, timeout;
  private MinimaxBoard b;
  private Random r;

  // NOTE: profiling
  private static int visited, series_found;

  // Fisher–Yates shuffle
  private void shuffle(MNKCell[] vec) {
    for (int i = vec.length - 1; i > 0; i--) {
      int index = r.nextInt(i + 1);
      // Simple swap
      MNKCell a = vec[index];
      vec[index] = vec[i];
      vec[i] = a;
    }
  }

  enum Action {
    MINIMIZE, MAXIMIZE
  }

  private static Action opposite(final Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  private static MNKCellState opponent(final MNKCellState p) {
    return p == MNKCellState.P1 ? MNKCellState.P2 : MNKCellState.P1;
  }

  private boolean should_halt() {
    // TODO: tweak values
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(95.0/100.0);
  }

  private double evaluate(final MNKBoard board, final int depth) {
    if(board.gameState() == MY_WIN)
      return RANK_CONSTANT / depth;
    else if(board.gameState() == ENEMY_WIN)
      return -(depth * RANK_CONSTANT);
    else if(board.gameState() == MNKGameState.DRAW)
      return 0;
    else throw new Error("invalid board state: game not ended");
  }

  private Pair<Double, MNKCell> minimax(MinimaxBoard board, Action action, int depth, double a, double b) {
    visited++;
    // handle the first move by placing ourselves at the center, which is the best postition for any mnk
    // if(board.getMarkedCells().length == 0)
    //   return new Pair<>(Double.MAX_VALUE, new MNKCell(N/2, M/2, ME));

    if(board.gameState() != MNKGameState.OPEN)
      // return the evaluation of the current board with the last marked cell
      // as the decision which has brough up to this game state
      return new Pair<>(evaluate(board, depth), board.getMarkedCells()[board.getMarkedCells().length-1]);

    if(should_halt())
      return new Pair<>(HALT, board.getFreeCells()[0]);

    Series.Result series = Series.findSeries(board, action == Action.MAXIMIZE ? ME : ENEMY);

    // instantly pick one-move win/loss prevention moves (k-1 seires, for both players)
    // as they are certainly the best for the game outcome
    if(series.k1.size() > 0) {
      // if we know we are acting on the root we can just return without
      // evaluating the value of the move, as this won't be used anywhere
      if(depth == 0)
        return new Pair<>(0d, series.k1.get(0));

      board.markCell(series.k1.get(0));
      Pair<Double, MNKCell> result = minimax(board, opposite(action), depth+1, a, b);
      board.unmarkCell();
      return new Pair<>(result.first, series.k1.get(0));
    }
    
    double best = action == Action.MAXIMIZE ? -Double.MAX_VALUE : Double.MAX_VALUE;
    MNKCell best_cell = null;
    // TODO: (?) shuffle(series.noseries);
    for(MNKCell c : series.kn) {
      board.markCell(c);
      Pair<Double, MNKCell> rank = minimax(board, opposite(action), depth+1, a, b);
      board.unmarkCell();
      if(rank.first == HALT)
        return rank;

      if((action == Action.MAXIMIZE && (double)rank.first > best) || (action == Action.MINIMIZE && rank.first < best)) {
        best = rank.first;
        best_cell = c;
      }
      if(action == Action.MAXIMIZE && best > a)
        a = best;
      else if(action == Action.MINIMIZE && best < b)
        b = best;

      if(b <= a)
        break;
    }
    return new Pair<>(best, best_cell);
  }

  public void initPlayer(int M, int N, int K, boolean first, int timeout_in_secs) {
    this.M = M;
    this.N = N;
    this.K = K;
    this.timeout = timeout_in_secs;

    MY_WIN   = first ? MNKGameState.WINP1 : MNKGameState.WINP2; 
    ENEMY_WIN = first ? MNKGameState.WINP2 : MNKGameState.WINP1;
    ME = first ? MNKCellState.P1 : MNKCellState.P2;
    ENEMY = first ? MNKCellState.P2 : MNKCellState.P1;
    b = new MinimaxBoard(M, N, K);
    r = new Random(0);
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    r.setSeed(start_time); // have a new pseudo-random generator each turn to improve randomness
    visited = series_found = 0;

    if(MC.length > 0)
      b.markCell(MC[MC.length-1]); // keep track of the opponent's marks

    try {
      Pair<Double, MNKCell> result = minimax(b, Action.MAXIMIZE, 0, -Double.MAX_VALUE, Double.MAX_VALUE);
      System.out.println(playerName() + "\t: visited " + visited + " nodes, ended with result: " + result);
      System.out.println(playerName() + "\t: found a total of " + series_found + " free cells in series (up to k-3)");

      if(FC.length != b.getFreeCells().length) {
        System.out.println("FATAL: minimax didn't clean the board");
        return FC[0];
      }

      b.markCell(result.second);
      return result.second;
    } catch(Exception e) {
      e.printStackTrace();
      return FC[0];
    }
  }

  public String playerName() {
    return "David Hilbert";
  }
}
