package gametree;
import java.util.Iterator;

class GameTreeIterator<P> implements Iterator<GameTree<P>> {
  GameTree<P> tree;
  int i = 0;
  GameTreeIterator(GameTree<P> tree) {
    this.tree = tree;
  }

  // chekcs if we can iterator forward
  public boolean hasNext() {
    return i >= tree.getChildCount();
  }
    
  // returns the next child
  public GameTree<P> next() {
    return tree.getChildren()[i++];
  }
}
