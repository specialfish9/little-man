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

import java.util.List;
import java.util.ArrayList;

import java.util.concurrent.Future;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeoutException;
import java.util.concurrent.Callable;


/**
 * Runs a game against two MNKPlayer classes and prints the game scores:
 * <ul>
 * <li> 3 if second player wins (and the first player is not interrupted)</li>
 * <li> 2 if first player wins or if the adversary is interrupted (illegal move or timeout interrupt)</li>
 * <li> 1 if the game ends in draw </li>
 * </ul>
 * <p>
 * Usage: MNKPlayerTester [OPTIONS] &lt;M&gt; &lt;N&gt; &lt;K&gt; &lt;MNKPlayer class name&gt; &lt;MNKPlayer class name&gt;<br/>
 * OPTIONS:<br>
 * &nbsp;&nbsp;-t &lt;timeout&gt; Timeout in seconds</br>
 * &nbsp;&nbsp;-r &lt;rounds&gt;  &nbsp;Number of rounds</br>
 * &nbsp;&nbsp;-v &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Verbose
 * </p>
 */
public class MNKPlayerTester {
	private static int     TIMEOUT = 10;
	private static int     ROUNDS  = 1;
	private static boolean VERBOSE = false;

	private static int     M;
	private static int     N;
	private static int     K;

	private static MNKBoard B;

	private static MNKPlayer[] Player = new MNKPlayer[2];


	/** Scoring system */
	private static int WINP1SCORE = 2;
	private static int WINP2SCORE = 3;
	private static int DRAWSCORE  = 1;
	private static int ERRSCORE   = 2;

	private enum GameState {
		WINP1, WINP2, DRAW, ERRP1, ERRP2;
	}
	

	private MNKPlayerTester() {
	}

	
	private static void initGame() {
		if(VERBOSE) System.out.println("Initializing " + M + "," + N + "," + K + " board");
		B = new MNKBoard(M,N,K);
		// Timed-out initializaton of the MNKPlayers
		for(int k = 0; k < 2; k++) {
			if(VERBOSE) if(VERBOSE) System.out.println("Initializing " + Player[k].playerName() + " as Player " + (k+1));
			final int i = k; // need to have a final variable here 
			final Runnable initPlayer = new Thread() {
				@Override 
				public void run() { 
					Player[i].initPlayer(B.M,B.N,B.K,i == 0,TIMEOUT);
				}
			};

			final ExecutorService executor = Executors.newSingleThreadExecutor();
			final Future future = executor.submit(initPlayer);
			executor.shutdown();
			try { 
				future.get(TIMEOUT, TimeUnit.SECONDS); 
			} 
			catch (TimeoutException e) {
				System.err.println("Error: " + Player[i].playerName() + " interrupted: initialization takes too much time");
				System.exit(1);
			}
			catch (Exception e) { 
				System.err.println(e);
				System.exit(1);		
			}
			if (!executor.isTerminated())
				executor.shutdownNow();
		}
		if(VERBOSE) System.out.println();
	}

	private static class StoppablePlayer implements Callable<MNKCell> {
		private final MNKPlayer P;
		private final MNKBoard  B;

		public StoppablePlayer(MNKPlayer P, MNKBoard B) {
			this.P = P;
			this.B = B;
		}

		public MNKCell call()  throws InterruptedException {
			return P.selectCell(B.getFreeCells(),B.getMarkedCells());
		}
	}

	private static GameState runGame() {
		while(B.gameState() == MNKGameState.OPEN) {
			int  curr = B.currentPlayer();
			final ExecutorService executor = Executors.newSingleThreadExecutor();
			final Future<MNKCell> task     = executor.submit(new StoppablePlayer(Player[curr],B));
			executor.shutdown(); // Makes the  ExecutorService stop accepting new tasks
			
			MNKCell c = null;
			
			try {
				c = task.get(TIMEOUT, TimeUnit.SECONDS);
			}
			catch(TimeoutException ex) {
				int n = 3; // Wait some more time to see if it stops
				System.err.println("Player " + (curr+1) + " (" + Player[curr].playerName() + ") interrupted due to timeout");
				while(!task.isDone() && n > 0) {
					System.err.println("Waiting for " + Player[curr].playerName() + " to stop ... (" + n + ")");
					try {Thread.sleep(TIMEOUT*1000);} catch(InterruptedException e) {}
					n--;
				}
				
				if(n == 0) {
					System.err.println("Player " + (curr+1) + " (" +Player[curr].playerName() + ") still running: game closed");
					System.exit(1);
				} else {
					System.err.println("Player " + (curr+1) + " (" + Player[curr].playerName() + ") eventually stopped: round closed");
					return curr == 0 ? GameState.ERRP1 : GameState.ERRP2; 
				}
			}
			catch (Exception ex) {
				int n = 3; // Wait some more time to see if it stops
				System.err.println("Player " + (curr+1) + " (" + Player[curr].playerName() + ") interrupted due to exception");
				System.err.println(" " + ex);
				while(!task.isDone() && n > 0) {
					System.err.println("Waiting for " + Player[curr].playerName() + " to stop ... (" + n + ")");
					try {Thread.sleep(TIMEOUT*1000);} catch(InterruptedException e) {}
					n--;
				}
				if(n == 0) {
					System.err.println("Player " + (curr+1) + " (" + Player[curr].playerName() + ") still running: game closed");
					System.exit(1);
				} else {
					System.err.println("Player " + (curr+1) + " (" + Player[curr].playerName() + ") eventually stopped: round closed");
					return curr == 0 ? GameState.ERRP1 : GameState.ERRP2;
				}
			}
			
			if (!executor.isTerminated())
				executor.shutdownNow();

			if(B.cellState(c.i,c.j) == MNKCellState.FREE) {
				if(VERBOSE) System.out.println("Player " + (curr+1) + " (" + Player[curr].playerName() + ") -> [" + c.i + "," + c.j + "]");
				B.markCell(c.i,c.j);
			} else {
				System.err.println("Player " + (curr+1) + " (" + Player[curr].playerName() + ")  selected an illegal move [" + c.i + "," + c.j + "]: round closed");
				return curr == 0 ? GameState.ERRP1 : GameState.ERRP2;
			}
		}
		
		return B.gameState() == MNKGameState.DRAW ? GameState.DRAW : (B.gameState() == MNKGameState.WINP1 ? GameState.WINP1 : GameState.WINP2);
	}

	private static void parseArgs(String args[]) {
		List<String> L = new ArrayList<String>(); 
		for (int i = 0; i < args.length; i++) {
			switch(args[i].charAt(0)) {
				case '-':
					char c = (args[i].length() != 2 ? 'x' : args[i].charAt(1));
					switch(c) {
						case 't': 
							if(args.length < i+2)
								throw new IllegalArgumentException("Expected parameter after " + args[i]);
							
							try {
								TIMEOUT = Integer.parseInt(args[++i]);
							} catch(NumberFormatException e) {
								throw new IllegalArgumentException("Illegal integer format for " + args[i-1] + " argument: " + args[i]);
							}
							break;
						case 'r':
							if(args.length < i+2)
								throw new IllegalArgumentException("Expected parameter after " + args[i]);	
							
							try {
								ROUNDS = Integer.parseInt(args[++i]);
							} catch(NumberFormatException e) {
								throw new IllegalArgumentException("Illegal integer format for " + args[i-1] + " argument: " + args[i]);
							}
							break;
						case 'v':
							VERBOSE = true;
							break;
						default: 
							throw new IllegalArgumentException("Illegal argument:  " + args[i]);
					}
					break;
				default:
				  L.add(args[i]);
			}
		}

		int n = L.size();
		if(n != 5)
			throw new IllegalArgumentException("Missing arguments:" + (n < 1 ? " <M>" : "") + (n < 2 ? " <N>" : "") +
			  (n < 3 ? " <K>" : "") + (n < 4 ? " <MNKPlayer class>" : "") + (n < 5 ? " <MNKPlayer class>" : "")); 

		try {
			M  = Integer.parseInt(L.get(0));
		}
		catch(NumberFormatException e) {
			throw new IllegalArgumentException("Illegal integer format for M argument: " + M);
		}
		try {
			N  = Integer.parseInt(L.get(1));
		}
		catch(NumberFormatException e) {
			throw new IllegalArgumentException("Illegal integer format for N argument: " + N);
		}
		try {
			K  = Integer.parseInt(L.get(2));
		}
		catch(NumberFormatException e) {
			throw new IllegalArgumentException("Illegal integer format for N argument: " + K);
		}

		if(M <= 0 || N <= 0 || K <= 0)
			throw new IllegalArgumentException("Arguments  M, N, K must be larger than 0");

		String[] P = {L.get(3),L.get(4)};
		for(int i = 0; i < 2; i++) {
			try {
				Player[i] = (MNKPlayer) Class.forName(P[i]).getDeclaredConstructor().newInstance();
			}
			catch(ClassNotFoundException e) {
				throw new IllegalArgumentException("Illegal argument: \'" + P[i] + "\' class not found");
      }
      catch(ClassCastException e) {
				throw new IllegalArgumentException("Illegal argument: \'" + P[i] + "\' class does not implement the MNKPlayer interface");
      }
      catch(NoSuchMethodException e) {
				throw new IllegalArgumentException("Illegal argument: \'" + P[i] + "\' class constructor needs to be empty");
      }
			catch(Exception e) {
				throw new IllegalArgumentException("Illegal argument: \'" + P[i] + "\' class (unexpected exception) " + e);
			}
    }
	}

	private static void printUsage() {
		System.err.println("Usage: MNKPlayerTester [OPTIONS] <M> <N> <K> <MNKPlayer class> <MNKPlayer class>");
		System.err.println("OPTIONS:");
		System.err.println("  -t <timeout>  Timeout in seconds. Default: " + TIMEOUT);
		System.err.println("  -r <rounds>   Number of rounds. Default: " + ROUNDS);
		System.err.println("  -v            Verbose. Default: " + VERBOSE);
	}

	public static void main(String[] args) {
		int P1SCORE = 0;
		int P2SCORE = 0;

		if(args.length == 0) {	
			printUsage();
			System.exit(0);
		}
		
		try {
			parseArgs(args);
		}
		catch(Exception e) {
			System.err.println(e);
			System.exit(1);	
		}
	
		if(VERBOSE) {
			System.out.println("Game type : " + M + "," + N + "," + K);
			System.out.println("Player1   : " + Player[0].playerName());
			System.out.println("Player2   : " + Player[1].playerName());
			System.out.println("Rounds    : " + ROUNDS);
			System.out.println("Timeout   : " + TIMEOUT + " secs\n\n");
		}

		for(int i = 1; i <= ROUNDS; i++) {
			if(VERBOSE) System.out.println("\n**** ROUND " + i + " ****");
			initGame();
			GameState state = runGame();

			switch(state) {
				case WINP1: P1SCORE += WINP1SCORE; break;
				case WINP2: P2SCORE += WINP2SCORE; break;
				case ERRP1: P2SCORE += ERRSCORE;   break;
				case ERRP2: P1SCORE += ERRSCORE;   break;
				case DRAW : P1SCORE += DRAWSCORE;
				            P2SCORE += DRAWSCORE;
				            break;
			}
			if(VERBOSE) {
				System.out.println("\nGame state    : " + state);
				System.out.println("Current score : " + Player[0].playerName() + " (" + P1SCORE + ") - " + Player[1].playerName() + " (" + P2SCORE + ")");
			}
		}
		if(VERBOSE) System.out.println("\n**** FINAL SCORE ****");
		System.out.println(Player[0].playerName() + " " + P1SCORE);
		System.out.println(Player[1].playerName() + " " + P2SCORE);	
	}
}
