package mnkgame.players;

import mnkgame.*;
import java.util.Random;



/*
 * AlphaBeta player with fixed depth minimax and montecarlo-like evaluation
 */
public class CharlesDarwin implements MNKPlayer {
  private final double RANK_CONSTANT = 10;
  private final double HALT = Double.MIN_VALUE+1;
  private MNKCellState ME;
  private MNKGameState MY_WIN, OTHER_WIN;

  private int M, N, K;
  private long start_time, timeout;
  private MinimaxBoard b;
  private int max_depth, simulations;

  private class Pair<A, B> {
    public A first;
    public B second;

    Pair(A first, B second) {
      this.first = first;
      this.second = second;
    }
  }

  enum Action {
    MINIMIZE, MAXIMIZE
  }

  private static Action opposite(Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  // simulates randomly a match starting from the given board n times
  // the ranking is done in such a way that it will work well with our alpha
  // beta pruning algorithm evaluation logic
  private double sim(MinimaxBoard board, int depth, final Random r, int n) {
    double wins = 0, losses = 0;
    while(n > 0) {
      int len = board.getFreeCells().length, played = 0;

      // make len random moves
      while(played < len && board.gameState() == MNKGameState.OPEN) {
        board.markCell(board.getFreeCells()[len-played-1]);
        played++;
      }

      // evaluate the random match
      if(board.gameState() == MY_WIN)
        wins++;
      else if(b.gameState() == MNKGameState.DRAW)
        wins += .5f; // we consider draws as wins, is that ok? (+= losses)
      else
        losses++;

      // cleanup all random moves
      while(played > 0) {
        b.unmarkCell();
        played--;
      }
      n--;
    }
    return wins > losses ? (RANK_CONSTANT / depth) * (wins/n) : -(depth * RANK_CONSTANT) * (losses/n);
  }

  private boolean should_halt() {
    // TODO: tweak values
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(99.0/100.0);
  }

  private Pair<Double, MNKCell> search(MinimaxBoard board, Action action, int depth, double a, double b, Random r) {
    // handle the first move by placin ourselves at the center, which is the best postition for any mnk
    if(board.getMarkedCells().length == 0)
      return new Pair<>(Double.MAX_VALUE, new MNKCell(N/2, M/2, ME));

    if(board.gameState() == MY_WIN)
      return new Pair<>(RANK_CONSTANT / depth, null);
    else if(board.gameState() == OTHER_WIN)
      return new Pair<>(-(depth * RANK_CONSTANT), null);
    else if(board.gameState() == MNKGameState.DRAW)
      return new Pair<>(0d, null);
    else if(depth >= max_depth)
      return new Pair<>(sim(board, depth, r, simulations), null);

    if(should_halt())
      return new Pair<>(HALT, board.getFreeCells()[0]);

    double best = action == Action.MAXIMIZE ? Double.MIN_VALUE : Double.MAX_VALUE;
    MNKCell best_cell = null;
    for(MNKCell c : board.getFreeCells()) {
      board.markCell(c);
      Pair<Double, MNKCell> rank = search(board, opposite(action), depth+1, a, b, r);
      board.unmarkCell();

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
    simulations = 3*K;
    max_depth = 3^K;
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    if(MC.length > 0)
      b.markCell(MC[MC.length-1]); // keep track of the opponent's marks

    try {
    Random r = new Random();
    Pair<Double, MNKCell> result = search(b, Action.MAXIMIZE, 0, Double.MIN_VALUE, Double.MAX_VALUE, r);
    b.markCell(result.second);
    return result.second;
    } catch(Exception e) {
      e.printStackTrace();
      return b.getFreeCells()[0];
    }
  }

  public String playerName() {
    return "Charles Darwin";
  }
}
