package mnkgame.players;

import java.util.Random;

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
  private Random r;

  // NOTE: profiling
  private int visited;

  private class Pair<A, B> {
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

  // Fisher–Yates shuffle
  private void shuffleArray(MNKCell[] vec)
  {
    for (int i = vec.length - 1; i > 0; i--)
    {
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

  // finds the first cell needed to copmlete a K-1 streak in any possible direction
  private MNKCell find_one_move_win(MinimaxBoard board, final MNKGameState win_state) {
    for(MNKCell c : board.getFreeCells()) {
      if(should_halt())
        return null;

      MNKGameState result = board.markCell(c);
      board.unmarkCell();
      if(result == win_state)
        return c;
    }
    return null;
  }

  private MNKCell pick_random_non_closing_cell(MinimaxBoard board, MNKCell previous) {
    for(MNKCell c : board.getFreeCells()) {
      // avoiding the should_halt check here, kinda superflous
      MNKGameState result = board.markCell(c);
      board.unmarkCell();
      if(result == MNKGameState.OPEN && (previous == null || previous.i != c.i || previous.j != c.j))
        return c;
    }
    return null;
  }

  // finds the first cell the enemy needs to copmlete a K-1 streak in any possible direction
  private MNKCell find_one_move_loss(MinimaxBoard board, final MNKGameState loss_state) {
    MNKCell random_cell = null;
    if(board.getFreeCells().length == 1 || (random_cell = pick_random_non_closing_cell(board, null)) == null)
      return null; // cannot check for enemy's next move when it doesn't exist

    board.markCell(random_cell);
    MNKCell c = find_one_move_win(board, loss_state);
    board.unmarkCell(); // remove the marked random_cell
    if(c != null)
      return c;

    // test the random_cell we selected at first. It may be a one-move loss cell
    // get a new random cell different from the previous and call it safe_cell
    if(board.markCell(pick_random_non_closing_cell(board, random_cell)) != MNKGameState.OPEN) {
      // random_cell puts us in a draw, ignore that
      board.unmarkCell();
      return null;
    }
    MNKGameState result = board.markCell(random_cell); // let the enemy take the random ane
    board.unmarkCell();
    board.unmarkCell();
    return result == loss_state ? random_cell : null;
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

    MNKCell omc = null; // one move win/loss cell
    if(board.getMarkedCells().length >= (K-1)*2 && // only check for one-win-moves
                                                   // if there have been placed
                                                   // enough cells to make one happen
      ((omc = find_one_move_win(board, action == Action.MAXIMIZE ? MY_WIN : ENEMY_WIN)) != null ||
      (omc = find_one_move_loss(board, action == Action.MAXIMIZE ? ENEMY_WIN : MY_WIN)) != null)) {
      // if we know we are acting on the root we can just return without
      // evaluating the value of the move, as this won't be used anywhere
      if(depth == 0)
        return new Pair<>(0d, omc);

      board.markCell(omc);
      Pair<Double, MNKCell> result = minimax(board, opposite(action), depth+1, a, b);
      board.unmarkCell();
      return new Pair<>(result.first, omc);
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
    r = new Random(0);
  }

  public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    start_time = System.currentTimeMillis();
    r.setSeed(start_time); // have a new pseudo-random generator each turn to improve randomness
    visited = 0;

    if(MC.length > 0)
      b.markCell(MC[MC.length-1]); // keep track of the opponent's marks

    try {
      Pair<Double, MNKCell> result = minimax(b, Action.MAXIMIZE, 0, -Double.MAX_VALUE, Double.MAX_VALUE);
      System.out.println(playerName() + "\t: visited " + visited + " nodes, ended with result: " + result);

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
    return "Charles Darwin";
  }
}
