package mnkgame.players;

import mnkgame.*;
import gametree.GameTree;

enum LeafResult {
  LOST,
  DRAW,
  VICTORY
}

/*
 * MiniMax player
 */
public class MinimaxPlayer implements MNKPlayer {
  int M, N, K;
  GameTree<MNKCell[]> tree;

  // generates a M * N board of cells, all set to be free
  private MNKCell[] freeState() {
    MNKCell[] grid = new MNKCell[M * N];
    for(int i = 0; i < M; i++)
      for(int j = 0; j < N; j++)
        grid[(i*N)+j] = new MNKCell(i, j);

    return grid;
  }

	public void initPlayer(int M, int N, int K, boolean first, int timeout_in_secs) {
    this.M = M;
    this.N = N;
    this.K = K;

    this.tree = new GameTree(freeState(), M*N);
    // run minimax and build the decision tree
    // ofc this is brutal dumb.
    // AT LEAST we ought to wait for the first move by the first player, then
    // start computing game trees
  }

	public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
    return null;
  }

	public String playerName() {
		return "MiniMax";
	}
}
