/*
 *  Copyright (C) 2021 Pietro Di Lena
 *  
 *  This file is part of the MNKGame v2.0 software developed for the
 *  students of the course "Algoritmi e Strutture di Dati" first 
 *  cycle degree/bachelor in Computer Science, University of Bologna
 *  A.Y. 2020-2021.
 *
 *  MNKGame is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This  is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this file.  If not, see <https://www.gnu.org/licenses/>.
 */

package mnkgame;

import java.lang.IndexOutOfBoundsException;
import java.lang.IllegalStateException;
import java.util.LinkedList;
import java.util.HashSet;

/**
 * Board for an (M,N,K)-game.
 * <p>
 * The MNKBoard class allows only alternates moves between two players. It mantains the ordered
 * list of moves and allows undoes.
 * 
 * </p> 
 */
public class MNKBoard {
	/**
   * Board rows
   */
	public final int M;
	/**
   * Board columns
   */
	public final int N;
  /**
   * Number of symbols to be aligned (horizontally, vertically, diagonally) for a win
   */
	public final int K;

	protected final MNKCellState[][]    B;
	protected final LinkedList<MNKCell> MC;  // Marked Cells
	protected final HashSet<MNKCell>    FC;  // Free Cells

	private final MNKCellState[] Player = {MNKCellState.P1,MNKCellState.P2};

	protected int          currentPlayer;   // currentPlayer plays next move

	protected MNKGameState gameState;       // game state
	
	/**
   * Create a board of size MxN and initialize the game parameters
   * 
   * @param M Board rows
	 * @param N Board columns
	 * @param K Number of symbols to be aligned (horizontally, vertically, diagonally) for a win
	 *
   * @throws IllegalArgumentException If M,N,K are smaller than  1
   */
	public MNKBoard(int M, int N, int K) throws IllegalArgumentException {
		if (M <= 0) throw new IllegalArgumentException("M cannot be smaller than 1");
		if (N <= 0) throw new IllegalArgumentException("N cannot be smaller than 1");
		if (K <= 0) throw new IllegalArgumentException("K cannot be smaller than 1");

		this.M  = M;
		this.N  = N;
		this.K  = K;

		B  = new MNKCellState[M][N];
		// Initial capacity large enough to assure load factor < 0.75
		FC = new HashSet<MNKCell>((int) Math.ceil((M*N) / 0.75)); 
		MC = new LinkedList<MNKCell>();

		reset();
	}

	/**
 	 * Resets the MNKBoard
	 */
	public void reset() {
		currentPlayer = 0;
		gameState     = MNKGameState.OPEN;
		initBoard();
		initFreeCellList();
		initMarkedCellList();
	}
	
	/**
	 * Returns the state of cell <code>i,j</code>
   *
   * @param i i-th row
   * @param j j-th column
   *
   * @return State of the <code>i,j</code> cell (FREE,P1,P2)
   * @throws IndexOutOfBoundsException If <code>i,j</code> are out of matrix bounds
	 */
	public MNKCellState cellState(int i, int j) throws IndexOutOfBoundsException {
		if(i < 0 || i >= M || j < 0 || j >= N)
			throw new IndexOutOfBoundsException("Indexes " + i + "," + j + " are out of matrix bounds");
		else
			return B[i][j];
	}

	/**
	 * Returns the current state of the game.
   *
   * @return MNKGameState enumeration constant (OPEN,WINP1,WINP2,DRAW)
   */
	public MNKGameState gameState() {
    return gameState;
  }

	/**
	 * Returns the id of the player allowed to play next move. 
	 *
	 * @return 0 (first player) or 1 (second player)
	 */
	public int currentPlayer() {
		return currentPlayer;	
	}
	
	/**
   * Marks the selected cell for the current player
   *
   * @param i i-th row
	 * @param j j-th column
   * 
	 * @return State of the game after the move
	 * 
	 * @throws IndexOutOfBoundsException If <code>i,j</code> are out of matrix bounds 
	 * @throws IllegalStateException If the game already ended or if <code>i,j</code> is not a free cell
   */
	public MNKGameState markCell(int i, int j) throws IndexOutOfBoundsException, IllegalStateException {
		if(gameState != MNKGameState.OPEN) {
			throw new IllegalStateException("Game ended!");
		} else if(i < 0 || i >= M || j < 0 || j >= N) {
			throw new IndexOutOfBoundsException("Indexes " + i +"," + j + " out of matrix bounds");
		} else if(B[i][j] != MNKCellState.FREE) {
			throw new IllegalStateException("Cell " + i +"," + j + " is not free");
		} else {
			MNKCell oldc = new MNKCell(i,j,B[i][j]);
			MNKCell newc = new MNKCell(i,j,Player[currentPlayer]);

			B[i][j] = Player[currentPlayer];

			FC.remove(oldc);
			MC.add(newc);
			
			currentPlayer = (currentPlayer + 1) % 2;

			if(isWinningCell(i,j))
				gameState =  B[i][j] == MNKCellState.P1 ? MNKGameState.WINP1 : MNKGameState.WINP2;
			else if(FC.isEmpty())
				gameState = MNKGameState.DRAW;
			
			return gameState;
		}
  }

	/**
   * Undoes last move
   *
   * @throws IllegalStateException If there is no move to undo
   */
	public void unmarkCell() throws IllegalStateException {
		if(MC.size() == 0) {
			throw new IllegalStateException("No move to undo");
		} else {
			MNKCell oldc = MC.removeLast();
			MNKCell newc = new MNKCell(oldc.i,oldc.j,MNKCellState.FREE);

			B[oldc.i][oldc.j] = MNKCellState.FREE;
			
			FC.add(newc);
			currentPlayer = (currentPlayer + 1) % 2;
			gameState     = MNKGameState.OPEN;
		}
	}

	/**
	 * Returns the marked cells list in array format.
	 * <p>This is the history of the game: the first move is in the
	 * array head, the last move in the array tail.</p>
	 * @return List of marked cells
	 */ 
	public MNKCell[] getMarkedCells() {
		return MC.toArray(new MNKCell[MC.size()]);
	}

	/**
	 * Returns the free cells list in array format.
	 * <p>There is not a predefined order for the free cells in the array</p>
	 * @return List of free cells
	 */
	public MNKCell[] getFreeCells() {
		return FC.toArray(new MNKCell[FC.size()]);
	}

	// Sets to free all board cells
	private void initBoard() {
		for(int i = 0; i < M; i++)
			for(int j = 0; j < N; j++)
				B[i][j] = MNKCellState.FREE;
	}

	// Rebuilds the free cells set 
	private void initFreeCellList() {
		this.FC.clear();
		for(int i = 0; i < M; i++)
			for(int j = 0; j < N; j++)
				this.FC.add(new MNKCell(i,j));
	}

	// Resets the marked cells list
	private void initMarkedCellList() {
		this.MC.clear();
	}

	// Check winning state from cell i, j
	private boolean isWinningCell(int i, int j) {
		MNKCellState s = B[i][j];
		int n;

		// Useless pedantic check
		if(s == MNKCellState.FREE) return false;

		// Horizontal check
		n = 1;
		for(int k = 1; j-k >= 0 && B[i][j-k] == s; k++) n++; // backward check
		for(int k = 1; j+k <  N && B[i][j+k] == s; k++) n++; // forward check   
		if(n >= K) return true;

		// Vertical check
		n = 1;
		for(int k = 1; i-k >= 0 && B[i-k][j] == s; k++) n++; // backward check
		for(int k = 1; i+k <  M && B[i+k][j] == s; k++) n++; // forward check
		if(n >= K) return true;
		

		// Diagonal check
		n = 1;
		for(int k = 1; i-k >= 0 && j-k >= 0 && B[i-k][j-k] == s; k++) n++; // backward check
		for(int k = 1; i+k <  M && j+k <  N && B[i+k][j+k] == s; k++) n++; // forward check
		if(n >= K) return true;

		// Anti-diagonal check
		n = 1;
		for(int k = 1; i-k >= 0 && j+k < N  && B[i-k][j+k] == s; k++) n++; // backward check
		for(int k = 1; i+k <  M && j-k >= 0 && B[i+k][j-k] == s; k++) n++; // backward check
		if(n >= K) return true;

		return false;
	}
}
