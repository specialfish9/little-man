package mnkgame.players;

import mnkgame.*;

/*
 * MiniMax player
 * Ideas:
 * - Generate the GameTree in place using the methods from MNKBoard, stopping
 *   when we find an enemy move which proves to be closing. (we always expect
 *   the enemy to be smart)
 * - 
 */
public class BertrandRussell implements MNKPlayer {
  private final double RANK_CONSTANT = 10;
  private final double HALT = Double.MIN_VALUE;
  private MNKCellState ME;
  private MNKGameState MY_WIN, OTHER_WIN;

  private int M, N, K;
  private long start_time, timeout;
  private MinimaxBoard b;

  // NOTE: profiling
  private int visited;

  private class Pair<A, B> {
    public A first;
    public B second;

    Pair(A first, B second) {
      this.first = first;
      this.second = second;
    }

    public String toString() {
      return "(" + first + "," + second + ")";
    }
  }

  enum Action {
    MINIMIZE, MAXIMIZE
  }

  private static Action opposite(Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  private boolean should_halt() {
    // TODO: tweak values
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(95.0/100.0);
  }

  private Pair<Double, MNKCell> minimax(MinimaxBoard board, Action action, int depth, double a, double b) {
    visited++;
    // MAYBE
    // handle the first move by placin ourselves at the center, which is the best postition for any mnk
    // if(board.getMarkedCells().length == 0)
    //  return new Pair<>(Double.MAX_VALUE, new MNKCell(N/2, M/2, ME));

    if(board.gameState() == MY_WIN)
      return new Pair<>(RANK_CONSTANT / depth, null);
    else if(board.gameState() == OTHER_WIN)
      return new Pair<>(-(depth * RANK_CONSTANT), null);
    else if(board.gameState() == MNKGameState.DRAW)
      return new Pair<>(0d, null);

    if(should_halt())
      return new Pair<>(HALT, board.getFreeCells()[0]);

    double best = action == Action.MAXIMIZE ? -Double.MAX_VALUE : Double.MAX_VALUE;
    MNKCell best_cell = null;
    for(MNKCell c : board.getFreeCells()) {
      board.markCell(c);
      Pair<Double, MNKCell> rank = minimax(board, opposite(action), depth+1, a, b);
      board.unmarkCell();
      if(rank.first == HALT)
        return rank;

      if(action == Action.MAXIMIZE && rank.first > best) {
        // during our turn take the best viable move
        best = a = rank.first;
        best_cell = c;
      } else if(action == Action.MINIMIZE && rank.first < best) {
        // during the opponent's turn we assume he takes the smartest move
        best = b = rank.first;
        best_cell = c;
      }

      if(b < a)
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
    OTHER_WIN = first ? MNKGameState.WINP2 : MNKGameState.WINP1;
    ME = first ? MNKCellState.P1 : MNKCellState.P2;
    b = new MinimaxBoard(M, N, K);
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    if(MC.length > 0)
      b.markCell(MC[MC.length-1]); // keep track of the opponent's marks

    visited = 0;
    Pair<Double, MNKCell> result = minimax(b, Action.MAXIMIZE, 0, -Double.MAX_VALUE, Double.MAX_VALUE);
    System.out.println(playerName() + "\t: visited " + visited + " nodes, ended with result: " + result);
    b.markCell(result.second);
    return result.second;
  }

  public String playerName() {
    return "Bertrand Russell";
  }
}
