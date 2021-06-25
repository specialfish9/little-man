package gametree;
import java.util.Iterator;

// A game-tree is a n-ary tree used to display the possibile plays in a
// trun-based game
//
// The implementation is kept generic, but it'll be mainly used to hold grid
// of MNKCells
//
// NOTE: we provide the maximum amount of children the tree will ever have, as
// we can know that at any state of the game. This allows us to use static arrays
// to lower the number of allocations.
public class GameTree<T> {
  private T state;
  private int maxNodes, nextNode = 0;
  private GameTree<T>[] children;

  public GameTree(T state, int maxNodes) {
    this.state = state;
    this.maxNodes = maxNodes;
    this.children = new GameTree[maxNodes];
  }
  public GameTree(T state) { this(state, 0); }

  public void setState(T state) {
    this.state = state;
  }

  public T getState() {
    return state;
  }

  public GameTree<T>[] getChildren() {
    return children;
  }

  public int getChildCount() {
    return nextNode;
  }

  public boolean isLeaf() {
    return maxNodes == 0 || children.length == 0;
  }

  public void append(T value) throws GameTreeFullException { append(value, 0); }
  public void append(T value, int maxNodes) throws GameTreeFullException {
    append(new GameTree(value, maxNodes));
  }

  public void append(GameTree node) throws GameTreeFullException {
    if(nextNode == maxNodes)
      throw new GameTreeFullException("Attempted to exceed the maximum amount" +
        " of children for the given GameTree");

    children[nextNode++] = node;
  }

  public Iterator<GameTree<T>> iterator() {
    return new GameTreeIterator<T>(this);
  }
}
