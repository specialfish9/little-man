package mnkgame.players;

import mnkgame.*;
import java.util.Random;

/*
 * Montecarlo player
 */
public class CharlesDarwin implements MNKPlayer {
  private final int HALT = Integer.MIN_VALUE+1;
  private final double C = 1.5;
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

  // simulates randomly a match starting from the given board
  private boolean sim(MinimaxBoard b, final Random r) {
    if(b.gameState() == MY_WIN)
      return true;
    else if(b.gameState() == OTHER_WIN)
      return false;
    else if(b.gameState() == MNKGameState.DRAW)
      return true; // TODO: se aggiorniamo il calcolo di x_i possiamo tenere meno conto delle draw

		MNKCell c = b.getFreeCells()[r.nextInt(b.getFreeCells().length)];
    b.markCell(c); // mark a random cell
    boolean res = sim(b, r);
    b.unmarkCell();
    return res;
  }

  // simulates `sims` times the game starting at the given MNKBoard.
  // Cost: \Theta(n*sims) where n is the number of empty cells
  private int get_xi(MinimaxBoard b, final Random r, int sims) {
    int wins = 0, i = sims;
    while(i-- > 0) {
      int free_cells = b.getFreeCells().length, j = free_cells; // n of cells to unmark later
      while(j-- > 0 && b.gameState() == MNKGameState.OPEN) {
        MNKCell c = b.getFreeCells()[r.nextInt(b.getFreeCells().length)];
        b.markCell(c);
      }

      if(b.gameState() == MY_WIN || b.gameState() == MNKGameState.DRAW)
        wins++;

      while(free_cells-- > 0)
        b.unmarkCell();
    }
    return (wins/sims)*100;
  }

  // upper confidence bound math
  private double uct(int xi, int N, double ni) {
    if(N == 0 || ni == 0)
      return Double.MAX_VALUE;
    else
      return xi + C*Math.sqrt(Math.log(N)/ni);
  }

  private final int SIMS = 5;
  public MNKCell montecarlo(MinimaxBoard b) {
    Random r = new Random();
    int free_cells = b.getFreeCells().length, N = 0;
    double[] vals = new double[free_cells];
    int max = 0;
    while(N < free_cells) {
      b.markCell(b.getFreeCells()[N]);
      vals[N] = uct(get_xi(b, r, SIMS), N, 0);
      b.unmarkCell();
      N++;
    }
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

    MNKCell result = montecarlo(b);
    b.markCell(result);
    return result;
  }

	public String playerName() {
		return "Charles Darwin";
	}
}
