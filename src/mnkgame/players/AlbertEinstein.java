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
public class AlbertEinstein implements MNKPlayer {
  private final int RANK_CONSTANT = 10;
  private final int HALT = Integer.MIN_VALUE+1;
  private MNKCellState ME;
  private MNKGameState MY_WIN, OTHER_WIN;

  private int M, N, K;
  private long start_time, timeout;
  private MinimaxBoard b;

  enum MatchResult {
    LOST,
    DRAW,
    VICTORY
  }

  // TODO: REMOVE
  private void printBoard(MinimaxBoard b) {
    for(int i = 0; i < M; i++) {
      for(int j = 0; j < N; j++) {
        System.out.print(b.cellState(i, j) + " ");
      }
      System.out.println("");
    }
  }

  private boolean should_halt() {
    // TODO: tweak values
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(99.0/100.0);
  }

  private int rank_board(MinimaxBoard b, boolean myTurn, int moves) {
    if((moves % 2 == 0 && should_halt()))
      return HALT;

    if(b.gameState() == MY_WIN)
      return RANK_CONSTANT / moves;
    else if(b.gameState() == OTHER_WIN)
      return -(moves * RANK_CONSTANT);
    else if(b.gameState() == MNKGameState.DRAW)
      return 0;

    int best = myTurn ? Integer.MIN_VALUE : Integer.MAX_VALUE;
    for(MNKCell c : b.getFreeCells()) {
      b.markCell(c);
      int rank = rank_board(b, !myTurn, moves+1);
      b.unmarkCell();
      if(myTurn)
        // during our turn take the best viable move
        if(rank > best)
          best = rank;
      else
        // during the opponent's turn we assume he takes the smartest move
        if(rank < best)
          best = rank;
    }
    return best;
  }

  // Alber(t) Einstein
  public MNKCell minimax(MinimaxBoard b) {
    /*
    let massimo -INF, cella_res = -1 
     
    for cellla : celle dispinibili
      let curr = calcola()
      if(curr > massimo) then
        massimo = curr
        cella_res = cella
      endif
    endfor
    return cella
*/
    if(b.getMarkedCells().length == 0)
      return new MNKCell(N/2, M/2, ME);

    int best_rank = Integer.MIN_VALUE;
    MNKCell best_cell = null; 
    System.out.println("Trovo ");
    for(MNKCell c : b.getFreeCells()) {
      b.markCell(c);
      int rank = rank_board(b, false, 1);
      b.unmarkCell();
      System.out.println(rank + " ");
      if(rank == HALT)
        return b.getFreeCells()[0]; // TODO: rethink

      if(rank > best_rank) {
        best_rank = rank;
        best_cell = c;
      }
    }
    System.out.print(" e prendo " + best_rank);
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
    System.out.println("I am the player" + Integer.toString(first ? 1 : 2));
  }

	public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    if(MC.length > 0)
			b.markCell(MC[MC.length-1]); // keep track of the opponent's marks
    printBoard(b);

    MNKCell result = minimax(b);
    System.out.println("marking " + result);
    b.markCell(result);
    return result;
  }

	public String playerName() {
		return "Albert Einstein";
	}
}
