import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;

public class EigenCFA {

//        #include "dirs.h"
//        #define COPY_STORE_ITERATION_NUMBER 50

  private final String PARAMS_PATH = "program-params.lst";
  private final String FUN_PATH = "Fun.mat";
  private final String ARG1_PATH = "Arg1.mat";
  private final String ARG2_PATH = "Arg2.mat";

  private final int SPARSE_FACTOR = 25;

  private final int COPY_STORE_ITERATION_NUMBER = 50;

  int STRIDE = 25;

  void populateMatrix(int[] matrix, int rows, int cols, String path) throws FileNotFoundException {
    Scanner scanner = new Scanner(new File(path)).useDelimiter("\\s|\\(|\\)");
    for (int i = 0; scanner.hasNext(); i++) {
      matrix[i] = scanner.nextInt();
    }
    scanner.close();
  }

  void reorder(int[] matrix, int count) {
    int i, j, initj;
    int[] temp = new int[count];
    int[] check = new int[count * 3];

    for (i = 0, initj = 0, j = 0; i < count; i++, j += STRIDE) {
      if (j >= count) {
        j = ++initj;
      }
      temp[i] = matrix[j];
      check[matrix[j]]++;
    }
    for (i = 0; i < count; i++) {
      if (check[temp[i]] < 1) {
        System.out.printf("Something really bad happened at %d\n", temp[i]);
      }
      matrix[i] = temp[i];
    }
  }

/*void populateStore (int* store, int rows, int cols, int lams, int WORD_SIZE)
{
  memset((void*)store, 0, rows * cols * sizeof(int));
  int i, count;
  for (i = rows*cols-1, count=0; count < lams; count++)
  {
	store[i] = (int)pow(2.0,WORD_SIZE-(i%WORD_SIZE)-1);
	if (i % WORD_SIZE == 0)
	  i -= (lams + 1);
	else
	  i--;
  }
  }*/

  void populateStore(int[][] store, int cols, int lams) {
    int vals = 3 * lams;

    for (int[] row : store) Arrays.fill(row, -1);

    for (int i = 0; i < lams; i++) {
      store[i][0] = 2;
      store[i][1] = i;
    }
    for (int i = lams; i < vals; i++) {
      store[i][0] = 1;
    }
  }


//  void displayMatrix(int*matrix, int rows, int cols) {
//    int i, j;
//    for (i = 0; i < rows; i++) {
//      for (j = 0; j < cols; j++) {
//        if (j % 10 == 0) {
//          printf("\n");
//        }
//        printf("%u ", * (matrix + i * cols + j));
//      }
//      printf("\n");
//    }
//  }


//  void displayStore(int*store, int rows, int cols) {
//    int i, j;
//    for (i = 0; i < rows; i++) {
//      int tail =*(store + i * cols + 0);
//      printf("%u - ", tail);
//      for (j = 1; j < tail; j++) {
//        printf("%d ", * (store + i * cols + j));
//      }
//      printf("\n");
//    }
//  }

  void reformatStore(int[][] store, int rows, int cols, File handle) throws FileNotFoundException {
    //  printf(" Lams : %d \n",lams);
    //  printf(" cvals: %d \n",cvals);
    PrintStream out = new PrintStream(handle);
    int i, j;
    out.printf("#hash(");
    for (i = 0; i < rows; i++) {
      int bindings = store[i][0] - 1;
      out.printf("(%d . (", i);
      for (j = 1; j <= bindings; j++) {
        out.printf("%d ", store[i][j]);
      }
      out.printf("))\n");
    }
    out.printf(")\n");
    out.close();
  }


  int find(int[] array, int start, int end, int searchTerm) {
    int i, elem;
    for (i = start; i < end && ((elem = array[i]) != -1); i++) {
      if (elem == searchTerm) {
        return -1;
      }
    }
    return i;
//  if (searchTerm == -1) {
// 	return -1;
//   }
//   int cell = searchTerm / WORDSIZE;
//   int bit = searchTerm % WORDSIZE;
//   return (array[cell] & (1 << (WORDSIZE-1 - bit)));
  }

  void updateStore(int[][] store_d, AtomicInteger nonZeroElements,
                   int[] callFun, int[] callArg1, int[] callArg2,
                   int lams, int calls, int cols, int initCall) {
    int i, j;
    int call;

    System.out.println("CALLS: " + calls);
    int old = nonZeroElements.get();

    for (call = 0; call < calls; call++) {

      int storeRowL = callFun[call];
      int storeRowL1 = callArg1[call];
      int storeRowL2 = callArg2[call];

      int lengthL = store_d[storeRowL][0];
      int lengthL1 = store_d[storeRowL1][0];
      int lengthL2 = store_d[storeRowL2][0];

      for (i = 1; i < lengthL; i++) {
        int row = store_d[storeRowL ][i];
        int v1Row = row + lams;
        int v2Row = row + 2 * lams;

        for (j = 1; j < lengthL1; j++) {
          int l1 = store_d[storeRowL1][j];
          int foundAt = find(store_d[v1Row], 1, cols, l1);
          if (foundAt != -1) {
            store_d[v1Row][store_d[v1Row][0]] = l1;
            store_d[v1Row][0]++;
            nonZeroElements.incrementAndGet();
          }
        }

        for (j = 1; j < lengthL2; j++) {
          int l2 = store_d[storeRowL2][j];
          int foundAt = find(store_d[v2Row], 1, cols, l2);
          if (foundAt != -1) {
            store_d[v2Row][store_d[v2Row][0]] = l2;
            store_d[v2Row][0]++;
            nonZeroElements.incrementAndGet();
          }
        }
      }
    }

    System.out.println(nonZeroElements.get() - old);
  }

  void eigenBitCFA(int[][] store,
                   int[] callFun, int[] callArg1, int[] callArg2,
                   int lams, int vals, int calls, int cols, int maxBindings,
                   int maxIterations) throws FileNotFoundException {
    float totalCompareTime = 0.0f, compareTime;
    float elapsedTime = 0.0f;
    float totalCopyTime = 0.0f, copyTime;

    // Global start and stop events for entire CFA
    long start, stop;
    // Events to measure time taken in comparing stores to find a fixed point
    long startCompare, stopCompare;
    // Events to measure the time to copy the updated store to the old store
    long startCopy, stopCopy;
    // Events to measure the time to copy the store from the device to the host
    long startCopyFromDevice, stopCopyFromDevice;

    int iterations;
    boolean iterate = true;
    AtomicInteger nonZeroElements_d = new AtomicInteger();
    int nonZeroElements = 0;
    int prevNonZeroElements = 0;
    System.out.printf("Starting analysis ... \n");

    // Start measuring total time to run the analysis
    start = System.currentTimeMillis();
    for (iterations = 0; iterate; iterations++) {

      updateStore(store, nonZeroElements_d,
          callFun, callArg1, callArg2, lams, calls, cols, 0);

      // Measure the time to compare the store
      startCompare = System.currentTimeMillis();

      nonZeroElements = nonZeroElements_d.get();
      iterate = (prevNonZeroElements != nonZeroElements);
      prevNonZeroElements = nonZeroElements;

      stopCompare = System.currentTimeMillis();
      compareTime = stopCompare - startCompare;
      totalCompareTime += compareTime;

      if (iterations % COPY_STORE_ITERATION_NUMBER == 0) {
        System.out.printf("iteration: %d (%d)\n", iterations, nonZeroElements);
      }
    }

    stop = System.currentTimeMillis();
    elapsedTime = stop - start;

    File handle = new File("scheme/tmp-stores/cpu_sparse_final");
    reformatStore(store, vals, cols, handle);

    System.out.printf("Analysis complete!\n");
    System.out.printf("\nTotal iterations: %d\n", iterations);
    System.out.printf("Total comparing time: %f\n", totalCompareTime / 1000);
    System.out.printf("Total copying time: %f\n", totalCopyTime / 1000);

    System.out.printf("\nTotal execution time: %f\n", elapsedTime / 1000);
    System.out.printf("-----------------------\n\n");
  }

  public void run(int stride) throws FileNotFoundException {
    STRIDE = stride;

    int maxIterations;

    Scanner scanner = new Scanner(new File(PARAMS_PATH));
    int lams = scanner.nextInt();
    int vars = scanner.nextInt();
    int calls = scanner.nextInt();
    int vals = scanner.nextInt();
    scanner.close();

    vals = 3 * lams;
    maxIterations = lams * vars;
    int maxBindings = (lams + SPARSE_FACTOR - 1) / SPARSE_FACTOR;
    int cols = maxBindings + 1;

    System.out.printf("Program parameters\n");
    System.out.printf("------------------\n");
    System.out.printf("lams: %d\nvars: %d\nvals: %d\ncalls: %d\n", lams, vars, vals, calls);
    System.out.printf("cols: %d\nmaxBindings: %d\n\n", cols, maxBindings);

    int[][] store = new int[vals][cols];
    int[] callFun = new int[calls];
    int[] callArg1 = new int[calls];
    int[] callArg2 = new int[calls];

    // Varying the strides to see what effect it has on execution time

    System.out.printf("STRIDE: %d\n", STRIDE);

    // Read in the initial store
    System.out.printf("Creating store (%d x %d) ... ", vals, cols);
    populateStore(store, cols, lams);
    System.out.printf("Populated store\n");
    //    displayStore(store, vals, cols);


    // Read in the FUN matrix
    System.out.printf("Reading CALLFUN (%d x %d) ... ", calls, 1);
    populateMatrix(callFun, calls, 2, FUN_PATH);
    System.out.printf("Populated FUN\n");
    reorder(callFun, calls);
    //	displayMatrix(callFun, 1, calls);

    // Read in the ARG1 matrix
    System.out.printf("Reading ARG1 (%d x %d) ... ", calls, 1);
    populateMatrix(callArg1, calls, 2, ARG1_PATH);
    reorder(callArg1, calls);
    System.out.printf("Populated ARG1\n");
    //    displayMatrix(callArg1, 1, calls);

    // Read in the ARG2 matrix
    System.out.printf("Reading ARG2 (%d x %d) ... ", calls, 1);
    populateMatrix(callArg2, calls, 2, ARG2_PATH);
    reorder(callArg2, calls);
    System.out.printf("Populated ARG2\n");
//    displayMatrix(callArg2, 1, calls);

    eigenBitCFA(store, callFun, callArg1, callArg2,
        lams, vals, calls, cols, maxBindings,
        maxIterations);
  }

}

