import java.io.File;
import java.io.FileNotFoundException;

public class Main {

  public static void main(String[] args) throws FileNotFoundException {
    System.out.println(new File("foo").getAbsolutePath());
    int stride = args.length > 0 ? Integer.parseInt(args[0]) : 25;
    EigenCFA cfa = new EigenCFA();
    cfa.run(stride);
  }

}
