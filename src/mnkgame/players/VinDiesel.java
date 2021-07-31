package mnkgame.players;

import java.util.Arrays;
import java.util.Optional;

import mnkgame.*;

/*
  TODO:
  1. Find a way to identify each move.
  2. If moves makes prune then add to the ads
  3. Chose the move acording to the ads
  4. Pray
  5. Test And Clear
*/

//https://core.ac.uk/download/pdf/225892254.pdf
public class VinDiesel implements MNKPlayer{
  private final double RANK_CONSTANT = 10;
  private final double HALT = Double.MIN_VALUE+1;
  private MNKCellState ME;
  private MNKGameState MY_WIN, OTHER_WIN;

  private final AdaptiveList maxADS = new AdaptiveList();
  private final AdaptiveList minADS = new AdaptiveList();


  private int M, N, K;
  private long start_time, timeout;
  private MinimaxBoard b;

  private class Pair<A, B> {
    public A first;
    public B second;

    Pair(A first, B second) {
      this.first = first;
      this.second = second;
    }
  }

  enum Action {
    MINIMIZE, MAXIMIZE
  }

  private static Action opposite(Action a) {
    return a == Action.MAXIMIZE ? Action.MINIMIZE : Action.MAXIMIZE;
  }

  private boolean should_halt() {
    // TODO: tweak values
    return (System.currentTimeMillis()-start_time)/1000.0 > timeout*(99.0/100.0);
  }

  private Pair<Double, MNKCell> minimax(MinimaxBoard board, Action action, int depth, double a, double b) {
    // handle the first move by placin ourselves at the center, which is the best postition for any mnk
    if(board.getMarkedCells().length == 0)
      return new Pair<>(Double.MAX_VALUE, new MNKCell(N/2, M/2, ME));

    if(board.gameState() == MY_WIN)
      return new Pair<>(RANK_CONSTANT / depth, null);
    else if(board.gameState() == OTHER_WIN)
      return new Pair<>(-(depth * RANK_CONSTANT), null);
    else if(board.gameState() == MNKGameState.DRAW)
      return new Pair<>(0d, null);

    if(should_halt())
      return new Pair<>(HALT, board.getFreeCells()[0]);

    double best = action == Action.MAXIMIZE ? Double.MIN_VALUE : Double.MAX_VALUE;
    MNKCell best_cell = null;

    if(action == Action.MAXIMIZE){
      for(int i = 0; i < maxADS.size; i ++){
        int code = maxADS.get(i);
        Optional<MNKCell> cell = getCellFromCode(code, board.getFreeCells());
        if(!cell.isPresent()) continue;
        // if is a valid move then do it
        board.markCell(cell.get());
        Pair<Double, MNKCell> rank = minimax(board, opposite(action), depth+1, a, b);
        board.unmarkCell();
        if (rank.first > best){
          best = a = rank.first;
          best_cell = cell.get();
        }
        if(b < a){
          maxADS.add(code);
          return new Pair<>(best, best_cell);
        }
      }
    } else {
      for(int i = 0; i < minADS.size; i ++){
        int code = minADS.get(i);
        Optional<MNKCell> cell = getCellFromCode(code, board.getFreeCells());
        if(!cell.isPresent()) continue;
        // if is a valid move then do it
        board.markCell(cell.get());
        Pair<Double, MNKCell> rank = minimax(board, opposite(action), depth+1, a, b);
        board.unmarkCell();
        if (rank.first > best){
          best = b = rank.first;
          best_cell = cell.get();
        }
        if(b < a){
          maxADS.add(code);
          return new Pair<>(best, best_cell);
        }
      }
    }
    
    // In case we didn't have any move in ads proceed with normal alpha beta

    for(MNKCell c : board.getFreeCells()) {
      board.markCell(c);
      Pair<Double, MNKCell> rank = minimax(board, opposite(action), depth+1, a, b);
      board.unmarkCell();

      if(action == Action.MAXIMIZE && rank.first > best) {
        // during our turn take the best viable move
        best = a = rank.first;
        best_cell = c;
      } else if(action == Action.MINIMIZE && rank.first < best) {
        // during the opponent's turn we assume he takes the smartest move
        best = b = rank.first;
        best_cell = c;
      }

      if(b < a){
        if(action == Action.MAXIMIZE)
          maxADS.add(getMoveCode(c));
        else
          minADS.add(getMoveCode(c));
        break;
      }
    }
    return new Pair<>(best, best_cell);
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

    Pair<Double, MNKCell> result = minimax(b, Action.MAXIMIZE, 0, Double.MIN_VALUE, Double.MAX_VALUE);
    b.markCell(result.second);
    return result.second;
  }

  public String playerName() {
    return "Vin Diesel";
  }

  private int getMoveCode(final MNKCell cell){
    return cell.i * 3 + cell.j;
  }

  private Optional<MNKCell> getCellFromCode(final int code, final MNKCell[] cells){
    return Arrays.stream(cells).takeWhile(cell -> getMoveCode(cell) == code).findFirst();
  }

  public class AdaptiveList{
    
    private Node head;

    public int size = 0;

    public void add(int move){
        if(head == null){
          head = new Node(move);
          size++;
          return;
        }
        Node ptr = head;
        if(ptr.value == move) return;
        while (ptr.next != null){
            
          if(ptr.next.value == move){
            // if the element is already in the list, swap it with its predecessor
            int tmp = ptr.next.value;
            ptr.next.value = ptr.value;
            ptr.value = tmp;
            return;
          }
          ptr = ptr.next;
        }
        // else append it at the end
        ptr.next = new Node(move);
        size ++;
    }

    @Override
    public String toString(){
      String s = "[";
      Node ptr = head;
      while(ptr != null){
        s += ptr.value + ", ";
        ptr = ptr.next;
      }
      s += "]";
      return s;
    }

    public int get(int index){
      int i = 0;
      Node ptr = head;
      while(i++ < index && ptr != null)
        ptr = ptr.next;
      if (ptr == null)
        throw new IndexOutOfBoundsException("Size is " + size + " and you tried to get " + index);
      return ptr.value;
    }

    private class Node {
      public int value;
      public Node next;

      public Node(int value){
        this.value = value;
      }
    }
  }
}