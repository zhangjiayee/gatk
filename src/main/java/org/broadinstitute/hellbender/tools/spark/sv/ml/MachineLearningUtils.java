package org.broadinstitute.hellbender.tools.spark.sv.ml;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * This is a base utils class that provides an API for packages that need to use classifiers. Specific classifiers can
 * be provided by extending GATKClassifier. Additional hyperparameter tuning strategies can be provided by extending
 * ClassifierTuner. These abstract base classes implement some common routines that should make adding new classifiers
 * or tuning strategies a little less tedious.
 */
public class MachineLearningUtils {


    /**
     * Form submatrix specified by selecting specified rows
     */
    public static RealMatrix sliceRows(final RealMatrix matrix, final int[] sliceRows) {
        final int[] allColumns = getRange(matrix.getColumnDimension());
        return matrix.getSubMatrix(sliceRows, allColumns);
    }


    /* here are multiple static functions that are generally useful for machine-learning tasks. Kept public to encourage re-use */

    /**
     * Return Integer array with elements {0, 1, 2, ..., numElements - 1}
     */
    public static Integer[] getRange(final Integer numElements) {
        if(numElements < 0) {
            throw new IllegalArgumentException("numElements must be >= 0");
        }
        final Integer[] range = new Integer[numElements];
        for(Integer i = 0; i < numElements; ++i) {
            range[i] = i;
        }
        return range;
    }

    /**
     * Return int array with elements {0, 1, 2, ..., numElements - 1}
     */
    public static int[] getRange(final int numElements) {
        if(numElements < 0) {
            throw new IllegalArgumentException("numElements must be >= 0");
        }
        return IntStream.range(0, numElements).toArray();
    }

    /**
     * Dereference array and return new array sampled from the array of supplied indices
     * @return arr[indices]
     */
    public static int[] slice(final int[] arr, final int[] indices) {
        final int[] sliced_arr = new int[indices.length];
        for(int i = 0; i < indices.length; ++i) {
            sliced_arr[i] = arr[indices[i]];
        }
        return sliced_arr;
    }

    /**
     * Dereference array and assign new values into the supplied indices: arr[indices] = newValues
     */
    public static void sliceAssign(final int[] arr, final int[] indices, final int[] newValues) {
        if(indices.length != newValues.length) {
            throw new IllegalArgumentException("length of indices does not match length of newValues");
        }
        for(int i = 0; i < indices.length; ++i) {
            arr[indices[i]] = newValues[i];
        }
    }

    /**
     * Find index to maximum value in array
     */
    public static int argmax(final double[] arr) {
        int bestIndex = -1;
        double bestValue = Double.MIN_VALUE;
        for(int index = 0; index < arr.length; ++index) {
            if(arr[index] > bestValue) {
                bestIndex = index;
                bestValue = arr[index];
            }
        }
        return bestIndex;
    }

    /**
     * Join to RealMatrices together so that columns of A come first, then columns of B.
     */
    public static RealMatrix concatenateColumns(final RealMatrix matrixA, final RealMatrix matrixB) {
        if(matrixA.getRowDimension() != matrixB.getRowDimension()) {
            throw new IllegalArgumentException("matrixA and matrixB do not have same number of rows.");
        }
        final RealMatrix matrixC = matrixA.createMatrix(
                matrixA.getRowDimension(),
                matrixA.getColumnDimension() + matrixB.getColumnDimension()
        );
        matrixC.setSubMatrix(matrixA.getData(), 0, 0);
        matrixC.setSubMatrix(matrixB.getData(), 0, matrixA.getColumnDimension());
        return matrixC;
    }

    /**
     * Return int array of indices that would sort supplied array, e.g. in pseudocode: arr[argsort(arr)] == sort(arr)
     */
    public static int[] argsort(final int[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparingInt(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
    }

    /**
     * Return int array of indices that would sort supplied array, e.g. in pseudocode: arr[argsort(arr)] == sort(arr)
     */
    public static <T extends Comparable<? super T>> int[] argsort(final T[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparing(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
    }

    /**
     * Return random permutation of elements {0, 1, ..., numElements - 1}, using Knuth shuffle algorithm.
     */
    public static int[] getRandomPermutation(final Random random, final int numElements) {
        // Knuth shuffle
        final int[] permutation = getRange(numElements);
        for(int i = numElements - 1; i > 0; --i) {
            final int swap_ind = random.nextInt(i);
            final int swap_val = permutation[swap_ind];
            permutation[swap_ind] = permutation[i];
            permutation[i] = swap_val;
        }
        return permutation;
    }

    /**
     * Geven predicted and correct labels, return proportion of predictions that are correct
     */
    public static double getPredictionAccuracy(final int[] predictedLabels, final int[] correctLabels) {
        int numCorrect = 0;
        for(int row = 0; row < correctLabels.length; ++row) {
            if(predictedLabels[row] == correctLabels[row]) {
                ++numCorrect;
            }
        }
        return numCorrect / (double)correctLabels.length;
    }

    /**
     * Convert stratify matrix to stratify array. Columns are binned and stratify array values are obtained by finding
     * unique rows.
     * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
     * @param stratifyMatrix matrix to convert to stratify array
     * @param numBins Number of bins into which to group values in each column
     * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
     *                                  the number of columns in consideration is reduced and the under-count rows are
     *                                  searched again for unique values.
     * @return stratifyArray, int array with length = number of rows in stratifyMatrix
     */
    public static int[] stratifyMatrixToStratifyArray(final RealMatrix stratifyMatrix, final int numBins,
                                                      final int minCountsPerStratifyValue) {
        return stratifyMatrixToStratifyArray(stratifyMatrix, numBins, minCountsPerStratifyValue, new HashSet<>());
    }

    /**
     * Convert stratify matrix to stratify array. Columns are binned and stratify array values are obtained by finding
     * unique rows. This version can specify some columns as categorical (unbinnable).
     * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
     * @param stratifyMatrix matrix to convert to stratify array
     * @param numBins Number of bins into which to group values in each column
     * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
     *                                  the number of columns in consideration is reduced and the under-count rows are
     *                                  searched again for unique values.
     * @param categoricalColumns Columns in categoricalColumns are not binned.
     * @return stratifyArray, int array with length = number of rows in stratifyMatrix
     */
    public static int[] stratifyMatrixToStratifyArray(final RealMatrix stratifyMatrix, final int numBins,
                                                      final int minCountsPerStratifyValue,
                                                      final Set<Integer> categoricalColumns) {
        // first form binned version of each column
        final int[][] binnedStratifyMatrix = new int[stratifyMatrix.getRowDimension()][stratifyMatrix.getColumnDimension()];
        for(int columnIndex = 0; columnIndex < stratifyMatrix.getColumnDimension(); ++columnIndex) {
            final int[] columnBins = categoricalColumns.contains(columnIndex) ?
                    getCategoryCodes(stratifyMatrix.getColumn(columnIndex))
                    : getBinnedColumn(stratifyMatrix.getColumn(columnIndex), numBins);
            for(int row = 0; row < stratifyMatrix.getRowDimension(); ++row) {
                binnedStratifyMatrix[row][columnIndex] = columnBins[row];
            }
        }

        // Get a set of unique rows
        Collection<List<Integer>> resultsToCheck = getUniqueRows(
                binnedStratifyMatrix,
                IntStream.range(0, stratifyMatrix.getRowDimension()).boxed().collect(Collectors.toList()),
                stratifyMatrix.getColumnDimension()
        );
        // While there are rows that have too few instances, decrease the number of columns under consideration (only
        // for those undersized rows)
        final Set<List<Integer>> uniqueResults = new HashSet<>();
        for(int useNumColumns = stratifyMatrix.getColumnDimension() - 1; useNumColumns > 0; --useNumColumns) {
            // Add undersized rows to a set to be reprocessed. Add properly sized rows to final uniqueResults set.
            final Set<List<Integer>> undersizedStratify = new HashSet<>();
            for(final List<Integer> uniqueResult : resultsToCheck) {
                if(uniqueResult.size() < minCountsPerStratifyValue) {
                    undersizedStratify.add(uniqueResult);
                } else {
                    uniqueResults.add(uniqueResult);
                }
            }
            if(undersizedStratify.isEmpty()) {
                break; // no problematic values, we're done
            } else if(useNumColumns == 1) {
                // can't consider fewer columns, lump remaining rows into the stratify value with the fewest members, to
                // make an "odd-ball" value.
                if(uniqueResults.isEmpty()) {
                    // handle edge case where every unique row has too few counts
                    uniqueResults.add(new ArrayList<>());
                }
                final List<Integer> smallest = uniqueResults.stream().min(Comparator.comparing(List::size))
                        .orElseThrow(NoSuchElementException::new);
                for(final List<Integer> uniqueResult : undersizedStratify) {
                    smallest.addAll(uniqueResult);
                }
                Collections.sort(smallest);
            } else {
                // find unique rows from subset that are too small, looking at one fewer column
                final List<Integer> tooSmall = undersizedStratify.stream().flatMap(Collection::stream).sorted()
                        .collect(Collectors.toList());
                resultsToCheck = getUniqueRows(binnedStratifyMatrix, tooSmall, useNumColumns);
            }
        }


        final int[] stratifyArray = new int[binnedStratifyMatrix.length];
        int stratifyValue = 0;
        for(final List<Integer> uniqueResult : uniqueResults) {
            for(final int index : uniqueResult) {
                stratifyArray[index] = stratifyValue;
            }
            ++stratifyValue;
        }
        return stratifyArray;
    }

    /**
     * Find unique values in column array, and return integer codes corresponding to those unique values
     * @return int array of length column.length, with values in range {0, 1, ..., numUniqueValues - 1}
     */
    private static int[] getCategoryCodes(final double[] column) {
        final int[] categoryCodes = new int[column.length];
        final Map<Double, Integer> codesMap = new HashMap<>();
        for(int i = 0; i < column.length; ++i) {
            final double value = column[i];
            final Integer code = codesMap.getOrDefault(value, null);
            if(code == null) {
                categoryCodes[i] = codesMap.size();
                codesMap.put(value, codesMap.size());
            } else {
                categoryCodes[i] = code;
            }
        }
        return categoryCodes;
    }

    /**
     * Bin column by
     * 1) choosing bins by sampling evenly in percentile space
     * 2) using bisection to map each value into its appropriate bin.
     * If NaN values are present, they are placed in highest bin.
     * @param column array of values to bin
     * @param numBins number of bins to distribute values into (e.g. number of unique integer values of output array)
     * @return int array with length column.length and values in {0, 1, ..., numBins-1}
     */
    private static int[] getBinnedColumn(final double[] column, final int numBins) {
        final int numUniqueColumnValues = (int)DoubleStream.of(column).distinct().count();
        if(numUniqueColumnValues <= numBins) {
            // Too much repetition to bin the data into the requested number of bins, just return unique values.
            return getCategoryCodes(column);
        }

        // Attempt to get requested number of percentiles, insisting that all percentiles are unique. If some values are
        // repeated, "percentiles" will not be evenly spaced, and it may not return exactly the requested number.
        // NOTE: the last percentile will not be used for binning (percentiles are "posts", bins are "fence") so request
        // one more percentile than bins
        final double[] percentiles = getUniquePercentiles(column, numBins + 1);

        // bin column to specified percentiles, dumping NaN into the last bin if it is present. Note that the percentiles
        // will be properly sized to account for the presence of NaN values.
        // NOTE: binarySearch is limited to ignore last percentile, because we don't want the max value in the array
        // being mapped past the maximum requested number of bins.
        return DoubleStream.of(column).mapToInt(
                v -> Double.isNaN(v) ? numBins - 1 : Arrays.binarySearch(percentiles, 0, percentiles.length, v)
        ).toArray();
    }

    /**
     * Return samples from Percentile object that are evenly-distributed in percentile space.
     * @param column array with values to sample from
     * @param numPercentiles Number of values to sample. If fewer than requested unique percentiles are found (if some
     *                       values repeat), sampling density will be increased to attempt to find the requested number.
     * @return array of length <= numPercentiles with values drawn from column that are approximately even in percentile
     *         space.
     */
    private static double[] getUniquePercentiles(final double[] column, final int numPercentiles) {
        final int numNaN = (int)Arrays.stream(column).filter(Double::isNaN).count();
        final int numSortablePercentiles = numNaN == 0 ? numPercentiles : numPercentiles - 1;
        final Percentile percentileEvaluator = new Percentile().withEstimationType(Percentile.EstimationType.R_1);
        percentileEvaluator.setData(column);
        double[] percentiles = percentileSpace(percentileEvaluator, numSortablePercentiles);
        if(percentiles.length == numSortablePercentiles) {
            return percentiles; // should be the case for typical non-repeating data
        }

        final int numSortable = column.length - numNaN;
        int numRequestHigh = numSortablePercentiles;
        while(percentiles.length < numSortablePercentiles && numRequestHigh < numSortable) {
            numRequestHigh = Math.min(2 * numRequestHigh, numSortable);
            percentiles = percentileSpace(percentileEvaluator, numRequestHigh);
        }
        if(percentiles.length == numSortablePercentiles) {
            return percentiles;
        }
        int numRequestLow = numSortablePercentiles;
        int range = numRequestHigh - numRequestLow;
        while(range > 1) {
            final int numRequest = numRequestLow + range / 2;
            percentiles = percentileSpace(percentileEvaluator, numRequest);
            if(percentiles.length < numSortablePercentiles) {
                numRequestLow = numRequest;
            } else if(percentiles.length > numSortablePercentiles){
                numRequestHigh = numRequest;
            } else {
                return percentiles;
            }
            range = numRequestHigh - numRequestLow;
        }

        // I'm not convinced that it's possible to get down to this point. It implies that you checked for one more
        // percentile but got 2 or more new unique values back. Just in case, insist on having *fewer* than requested,
        // since other approximate cases always yield that outcome.
        percentiles = percentileSpace(percentileEvaluator, numRequestLow);

        return percentiles;
    }

    /**
     * Return samples from Percentile object that are evenly-distributed in percentile space.
     * @param percentileEvaluator Percentile object that previously had setData() called to assign data.
     * @param numPercentiles >= 2
     */
    private static final double[] percentileSpace(final Percentile percentileEvaluator, final int numPercentiles) {
        if(numPercentiles < 2) {
            throw new IllegalArgumentException("numPercentiles must be >= 2");
        }
        final double low = 50.0 / percentileEvaluator.getData().length;
        final double high = 100.0 - low;
        final double coef = (high - low) / (numPercentiles - 1);
        return IntStream.range(0, numPercentiles).mapToDouble(i -> percentileEvaluator.evaluate(low + i * coef))
                .distinct().toArray();
    }


    /**
     * Find unique rows, and return set with indices to rows in original data
     * @param binnedStratifyMatrix  int[][] representation of data matrix
     * @param checkRows       only look for unique rows indexed by checkRows
     * @param useNumColumns   when checking for uniqueness, only consider columns in {0, 1, , ..., numUseColumns - 1}
     * @return set where each element is a list of row indices (from checkRows) so that for each row in the list,
     *         binnedStratifyMatrix[row] is identical.
     */
    private static Collection<List<Integer>> getUniqueRows(final int[][] binnedStratifyMatrix,
                                                           final List<Integer> checkRows,
                                                           final int useNumColumns) {
        // Use columns 0 to numUseColumns - 1, convert to string, and use as a key to collect indices into lists of unique rows
        return checkRows.stream().collect(Collectors.groupingBy(
                rowIndex -> Arrays.toString(Arrays.copyOfRange(binnedStratifyMatrix[rowIndex], 0, useNumColumns))
        )).values();
    }
}
