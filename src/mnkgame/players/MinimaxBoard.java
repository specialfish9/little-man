package mnkgame.players;

import mnkgame.*;

public class MinimaxBoard extends MNKBoard implements Cloneable {
  public MinimaxBoard(int M, int N, int K) {
    super(M, N, K);
  }

  public MNKGameState markCell(MNKCell c) throws IndexOutOfBoundsException, IllegalStateException {
    return markCell(c.i, c.j);
  }

  public Object clone(){
    try {
      return super.clone();
    } catch (CloneNotSupportedException e) {
      e.printStackTrace();
      return null;
    }
  }
}
