package mnkgame.players;

public abstract class LoggedPlayer {

  protected int _visited;
  protected int _cutoff;

  public String GetResults() {
    return "VISITED: " + _visited + " CUTOFF: " + _cutoff;
  }
}
