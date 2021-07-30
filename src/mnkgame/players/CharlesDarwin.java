package mnkgame.players;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Collections;
import mnkgame.*;

/*
 * Minimax + Alpha Beta + random selection of next move
 */
public class CharlesDarwin implements MNKPlayer {
  private final double RANK_CONSTANT = 10;
  private final double HALT = Double.MIN_VALUE+1;
  private MNKCellState ME, ENEMY;
  private MNKGameState MY_WIN, OTHER_WIN;

  private int M, N, K;
  private long start_time, timeout;
  private MinimaxBoard b;

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

  private boolean should_halt() {
    // TODO: tweak values
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(99.0/100.0);
  }

  // finds the cells needed to copmlete a K-1 streak in any possible direction
  // (vertical, horizontal, diagonal, counterdiagonal) for the ai
  private List<MNKCell> find_one_move_win_cells(MinimaxBoard board) {
    List<MNKCell> cells = new ArrayList<>();
    for(MNKCell c : board.getFreeCells()) {
      if(should_halt())
        return cells;
      else if(board.markCell(c) == MY_WIN)
        cells.add(c);
      board.unmarkCell();
    }
    return cells;
  }

  private Pair<Double, MNKCell> minimax(MinimaxBoard board, Action action, int depth, double a, double b) {
    // handle the first move by placin ourselves at the center, which is the best postition for any mnk
    if(board.getMarkedCells().length == 0)
      return new Pair<>(Double.MAX_VALUE, new MNKCell(N/2, M/2, ME));

    if(board.gameState() == MY_WIN)
      return new Pair<>(RANK_CONSTANT / depth, null);
    else if(board.gameState() == OTHER_WIN)
      return new Pair<>(-(depth * RANK_CONSTANT), null);
    else if(board.gameState() == MNKGameState.DRAW)
      return new Pair<>(0d, null);

    if(should_halt())
      return new Pair<>(HALT, board.getFreeCells()[0]);

    double best = action == Action.MAXIMIZE ? Double.MIN_VALUE : Double.MAX_VALUE;
    MNKCell best_cell = null;
    List<MNKCell> free_cells = new ArrayList<>(Arrays.asList(board.getFreeCells()));
    List<MNKCell> one_move_win_cells = find_one_move_win_cells(board);
    Collections.shuffle(free_cells); // \Theta(n) i guess, both ignorable still
    free_cells.removeAll(one_move_win_cells); // remove any overlapping cells
    while(true) {
      // select a cell, first from the one_move_win_cells collection,
      // later from the free_cells collection
      MNKCell c;
      if(one_move_win_cells.size() > 0)
        c = one_move_win_cells.remove(0);
      else if(free_cells.size() > 0)
        c = free_cells.remove(0);
      else break;

      board.markCell(c);
      Pair<Double, MNKCell> rank = minimax(board, opposite(action), depth+1, a, b);
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

      if(b < a) {
        System.out.println("cutoff, " + best_cell);
        break;
      }
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
    ME = first ? MNKCellState.P2 : MNKCellState.P1;
    b = new MinimaxBoard(M, N, K);
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    if(MC.length > 0)
      b.markCell(MC[MC.length-1]); // keep track of the opponent's marks

    Pair<Double, MNKCell> result = minimax(b, Action.MAXIMIZE, 0, Double.MIN_VALUE, Double.MAX_VALUE);
    b.markCell(result.second);
    return result.second;
  }

  public String playerName() {
    return "Charles Darwin";
  }
}
