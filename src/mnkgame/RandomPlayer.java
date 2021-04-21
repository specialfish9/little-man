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

import java.util.Random;

/**
 * Totally random software player.
 */
public class RandomPlayer  implements MNKPlayer {
	private Random rand;
	private int TIMEOUT;

	/**
   * Default empty constructor
   */
	public RandomPlayer() {
	}

	public void initPlayer(int M, int N, int K, boolean first, int timeout_in_secs) {
		// New random seed for each game
		rand    = new Random(System.currentTimeMillis()); 
		// Save the timeout for testing purposes
		TIMEOUT = timeout_in_secs;
		
		// Uncomment to chech the initialization timeout
		/* 
		try {
			Thread.sleep(1000*2*TIMEOUT);
		} 
		catch(Exception e) {
		}
		*/
	}

	/**
   * Selects a random cell in <code>FC</code>
   */
	public MNKCell selectCell(MNKCell[] FC, MNKCell[] MC) {
		// Uncomment to chech the move timeout
		/* 
		try {
			Thread.sleep(1000*2*TIMEOUT);
		} 
		catch(Exception e) {
		}
		*/
		return FC[rand.nextInt(FC.length)];
	}

	public String playerName() {
		return "R4nd0m";
	}
}
