package mnkgame.players;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Collections;
import java.util.LinkedHashSet;

import mnkgame.*;

/*
 * Minimax + Alpha Beta + one-win-euristics (next-move wind/loss interception)
 */
public class CharlesDarwin implements MNKPlayer {
  private final double RANK_CONSTANT = 10;
  private final double HALT = Double.MIN_VALUE;
  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, ENEMY_WIN;

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
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(99.0/100.0);
  }

  private double evaluate(MNKBoard board, int depth) {
    if(board.gameState() == MY_WIN)
      return RANK_CONSTANT / depth;
    else if(board.gameState() == ENEMY_WIN)
      return -(depth * RANK_CONSTANT);
    else if(board.gameState() == MNKGameState.DRAW)
      return 0;
    else throw new Error("invalid board state: game not ended");
  }

  // finds the first cell needed to copmlete a K-1 streak in any possible direction
  private MNKCell find_one_move_win(MinimaxBoard board, MNKGameState win_state) {
    for(MNKCell c : board.getFreeCells()) {
      if(should_halt())
        return null;
      else if(board.markCell(c) == win_state) {
        board.unmarkCell();
        return c;
      } else
        board.unmarkCell();
    }
    return null;
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

    MNKCell omwc = null; // one move win cell
    if((omwc = find_one_move_win(board, action == Action.MAXIMIZE ? MY_WIN : ENEMY_WIN)) != null) {
      board.markCell(omwc);
      double value = evaluate(board, depth+1);
      board.unmarkCell();
      return new Pair<>(value, omwc);
    }

    double best = action == Action.MAXIMIZE ? -Double.MAX_VALUE : Double.MAX_VALUE;
    MNKCell best_cell = null;
    for(MNKCell c : board.getFreeCells()) {
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
    return "Charles Darwin";
  }
}
