import java.util.Scanner;

public class Main {
    static final int CROSS = -1;
    static final int EMPTY = 0;
    static final int CIRCLE = 1;

    static int minmax(int[][] board, int m, int n, int depth, boolean isMax){
        if(isGameEnded(board, m, n)){
            if(hasWin(board, m, n, CROSS))
                return 10 - depth;
            else return -10 + depth;
        }
        int res = isMax ? Integer.MIN_VALUE : Integer.MAX_VALUE;
        for(int i = 0; i < n; i ++){
            for(int j = 0; j < m; j++){
                if(board[i][j] != EMPTY)
                    continue;
                board[i][j] = CROSS;
                int call = minmax(board, m, n, depth + 1, !isMax);
                board[i][j] = EMPTY;
                res = isMax ? Math.max(res, call) : Math.min(res, call);
            }
        }
        return res;
    }

    static int[][] makeAiMove(int[][] board, int n, int m){
        int bestI = -1, bestJ = -1;
        int res = Integer.MIN_VALUE;
        for(int i = 0; i < n; i ++){
            for(int j = 0; j < m; j++){
                if(board[i][j] != EMPTY)
                    continue;
                board[i][j] = CROSS;
                int call = minmax(board, m, n,  1, false);
                if(call > res){
                    res = call;
                    bestI = i;
                    bestJ = j;
                }
                board[i][j] = EMPTY;
            }
        }
        board[bestI][bestJ] = CROSS;
        return board;
    }

    static boolean isGameEnded(int[][] board, int m, int n){
        for(int i = 0; i < m; i++)
            for(int j = 0; j < n; j++)
                if(board[i][j] == EMPTY) return false; 
        return true;
    }

    static boolean hasWin(int[][] board, int m, int n, int player){
        for(int i = 0; i < n; i ++){
            boolean eq = true;
            for(int j = 0; j < m; j++)
                eq = eq && board[i][j] == player;
            if(eq)
                return true;
        }
        for(int i = 0; i < m; i ++){
            boolean eq = true;
            for(int j = 0; j < n; j++)
                eq = eq && board[j][i] == player;
            if(eq)
                return true;
        }
        return false;
    }

    static void printBoard(int[][] board, int m, int n){
        System.out.println("----------------------------");
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++) {
                    System.out.print(board[i][j] == CROSS ? "X" : board[i][j] == CIRCLE ? "O" : ".");
            }
            System.out.println();
        }
        System.out.println("----------------------------");
    }

    public static void main(String args[]){
        int[][] board = new int[3][3];
        for(int i = 0; i < 3; i ++)
            for(int j = 0; j < 3; j++)
                board[i][j] = EMPTY;
        
        
        Scanner s = new Scanner(System.in);
        boolean ai = false;    
        while(!isGameEnded(board, 3, 3)){
            if(ai){
                board = makeAiMove(board, 3, 3);
                ai = false;
            } else {
                int x = s.nextInt();
                int y = s.nextInt();
                board[x][y] = CIRCLE;
                ai = true;
            }
            printBoard(board, 3,3);
        }
        s.close();

    }
}