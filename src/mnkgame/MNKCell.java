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
 * Describes the state of a cell in the <code>MNKBoard</code>.
 * 
 * @see MNKBoard MNKBoard
 */
public class MNKCell {
	/**
   * Cell row index
   */
	public final int  i;
	/**
   * Cell column index
   */
	public final int  j;
	/**
   * Cell state
   */
	public final MNKCellState state;
	

	/**
   * Allocates a cell 
	 * 
   * @param i cell row index
   * @param j cell column index
   * @param state cell state 
   */
	public MNKCell(int i, int j, MNKCellState state) {
		this.i     = i;
		this.j     = j;
		this.state = state;
	}

	/**
	 * Allocates a free cell
	 *
	 * @param i cell row index
	 * @param j cell column index
	 * 
	 */
	public MNKCell(int i, int j) {
		this.i     = i;
		this.j     = j;
		this.state = MNKCellState.FREE;
	}

	@Override
	public boolean equals(Object o) {
		if(o == null)               return false;
		if(o == this)               return true; 	
		if(!(o instanceof MNKCell)) return false;

		MNKCell c = (MNKCell) o;

		return this.i == c.i && this.j == c.j && this.state == c.state;
	}

	@Override
	public int hashCode() {
		return this.toString().hashCode();
	}

	@Override
	public String toString() {
		return "[" + this.i + "," + this.j + "] -> " + this.state;
	}
}
