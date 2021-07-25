package mnkgame.players;

import mnkgame.*;

/*
 * AlphaBeta player
 */
public class BertrandRussell implements MNKPlayer {
  private final int RANK_CONSTANT = 10;
  private final int HALT = Integer.MIN_VALUE+1;
  private MNKCellState ME;
  private MNKGameState MY_WIN, OTHER_WIN;

  private int M, N, K;
  private long start_time, timeout;
  private MinimaxBoard b;

  enum Action {
    MINIMIZE, MAXIMIZE
  }

  private boolean should_halt() {
    // TODO: tweak values
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(99.0/100.0);
  }

  private int rank_board(MinimaxBoard board, Action action, int moves, int a, int b) {
    if((moves % 2 == 0 && should_halt()))
      return HALT; // TODO: change logic, use another way to stop the execution

    if(board.gameState() == MY_WIN)
      return RANK_CONSTANT / moves;
    else if(board.gameState() == OTHER_WIN)
      return -(moves * RANK_CONSTANT);
    else if(board.gameState() == MNKGameState.DRAW)
      return 0;

    int best = action == Action.MAXIMIZE ? Integer.MIN_VALUE : Integer.MAX_VALUE;
    for(MNKCell c : board.getFreeCells()) {
      board.markCell(c);
      int rank = rank_board(board, action == Action.MAXIMIZE ? Action.MINIMIZE
        : Action.MAXIMIZE, moves+1, a, b);
      board.unmarkCell();

      if(action == Action.MAXIMIZE && rank > best) {
        // during our turn take the best viable move
        best = a = rank;
      } else if(rank < best){
        // during the opponent's turn we assume he takes the smartest move
        best = b = rank;
      }

      if(b < a)
        break;
    }
    return best;
  }

  public MNKCell minimax(MinimaxBoard b) {
    if(b.getMarkedCells().length == 0)
      return new MNKCell(N/2, M/2, ME);

    int best_rank = Integer.MIN_VALUE;
    MNKCell best_cell = null; 
    for(MNKCell c : b.getFreeCells()) {
      b.markCell(c);
      int rank = rank_board(b, Action.MINIMIZE, 1, Integer.MIN_VALUE, Integer.MAX_VALUE);
      b.unmarkCell();
      if(rank == HALT)
        return b.getFreeCells()[0]; // TODO: rethink

      if(rank > best_rank) {
        best_rank = rank;
        best_cell = c;
      }
    }
    return best_cell;
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

    MNKCell result = minimax(b);
    b.markCell(result);
    return result;
  }

	public String playerName() {
		return "Bertrand Russell";
	}
}
