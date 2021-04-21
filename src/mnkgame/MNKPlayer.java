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

/**
 * Interface for a (M,N,K)-game software player.
 * <p>
 * The implementing classes need to provide a constructor that takes no arguments. The MNKPlayer
 * is initialized through the <code>initPlayer</code> method.
 * </p>
 */
public interface MNKPlayer {
	/**
   * Initialize the (M,N,K) Player
   *
   * @param M Board rows
   * @param N Board columns
   * @param K Number of symbols to be aligned (horizontally, vertically, diagonally) for a win
   * @param first True if it is the first player, False otherwise
	 * @param timeout_in_secs Maximum amount of time (in seconds) for selectCell 
   */
	public void initPlayer(int M, int N, int K, boolean first, int timeout_in_secs);
	
	/**
	 * Select a position among those listed in the <code>FC</code> array
	 *
	 * @param FC Free Cells: array of free cells
	 * @param MC Marked Cells: array of already marked cells, ordered with respect
   * to the game moves (first move is in the first position, etc)
   *
   * @return an element of <code>FC</code>
	 */
	public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC);	

	/**
   * Returns the player name
   *
	 * @return string 
   */
	public String playerName();
}
