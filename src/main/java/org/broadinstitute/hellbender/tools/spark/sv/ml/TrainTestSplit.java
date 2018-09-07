package org.broadinstitute.hellbender.tools.spark.sv.ml;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;

/**
 * Class that holds partitions of data sets (by row). trainRows and testRows are row indices into original data set,
 * appropriate for passing to sliceRows()
 */
public class TrainTestSplit {
    private final int[] trainRows;
    private final int[] testRows;

    TrainTestSplit(final int[] trainRows, final int[] testRows) {
        this.trainRows = trainRows;
        this.testRows = testRows;
    }


    /**
     * Static method to split data into two sets: training set and testing set.
     * @param trainingFraction proportion in [0, 1] of original data that should be placed in training set.
     * @param numRows          number of rows of original data
     * @param random           random number generator
     * @param stratify         Vector of integers listing stratification class. Data will be split so that each
     *                         stratify value has balanced number of instances in each set (e.g. proportion of
     *                         rows with stratify == 3 in training set will be approximately trainingFraction). If
     *                         null, split randomly. It is highly desirable to stratify at least based on class label.
     * @return TrainTestSplit with appropriate assigned indices.
     */
    public static TrainTestSplit getTrainTestSplit(final double trainingFraction, final int numRows,
                                                   final Random random, final int[] stratify) {
        if(stratify != null && stratify.length != numRows) {
            throw new IllegalArgumentException(
                    "stratify.length (" + stratify.length + ") != numRows (" + numRows + ")"
            );
        }

        final long numTrain = Math.round(numRows * trainingFraction);
        final long numTest = numRows - numTrain;
        if(numTrain <= 0) {
            if(trainingFraction < 0) {
                throw new IllegalArgumentException(
                        "trainingFraction (" + trainingFraction + ") must be in range [0, 1]"
                );
            }
            return new TrainTestSplit(new int[0], MachineLearningUtils.getRange(numRows));
        } else if(numTest <= 0) {
            if(trainingFraction > 1) {
                throw new IllegalArgumentException(
                        "trainingFraction (" + trainingFraction + ") must be in range[0, 1]"
                );
            }
            return new TrainTestSplit(MachineLearningUtils.getRange(numRows), new int[0]);
        }
        final int[] split_index_ordering = stratify == null ?
                MachineLearningUtils.getRandomPermutation(random, numRows)
                : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);
        final int[] trainRows = new int[(int)numTrain];
        final int[] testRows = new int[(int)numTest];
        int nextTrainInd = 1;
        int nextTestInd = 1;
        for(final int split_index : split_index_ordering) {
            if(numTrain * nextTestInd >= numTest * nextTrainInd) {
                // training set gets next index
                trainRows[nextTrainInd - 1] = split_index;
                ++nextTrainInd;
            } else {
                testRows[nextTestInd - 1] = split_index;
                ++nextTestInd;
            }
        }
        Arrays.sort(trainRows);
        Arrays.sort(testRows);
        return new TrainTestSplit(trainRows, testRows);
    }

    /**
     * Static method to get iterator over partitions of data, for implementing k-fold cross-validation.
     * @param numCrossvalidationFolds >= 2
     * @param numRows          number of rows in original data set
     * @param random           random number generator
     * @param stratify         Vector of integers listing stratification class. Cross-validation will split data so
     *                         that each stratify value has balanced number of instances in each fold. If null,
     *                         split randomly. It is highly desirable to stratify at least based on class label.
     * @return
     */
    public static Iterator<TrainTestSplit> getCrossvalidationSplits(final int numCrossvalidationFolds, final int numRows,
                                                                    final Random random, final int[] stratify) {
        if(numCrossvalidationFolds < 2) {
            throw new IllegalArgumentException("numCrossvalidationFolds (" + numCrossvalidationFolds + ") must be >= 2");
        }
        if(stratify != null && stratify.length != numRows) {
            throw new IllegalArgumentException(
                    "stratify.length (" + stratify.length + ") != numRows (" + numRows + ")"
            );
        }
        final int[] split_index_ordering = stratify == null ?
                MachineLearningUtils.getRandomPermutation(random, numRows)
                : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);

        return new FoldSplitIterator(split_index_ordering, numCrossvalidationFolds);
    }

    /**
     * Get ordering of data so that evenly distributing rows according to this order causes sampling that is
     * a) balanced in each stratify value
     * b) randomly-ordered within each stratify value
     */
    private static int[] getStratfiedIndexOrdering(final Random random, final int[] stratify) {
        /*
        logical (but memory inefficient) process
        1. make a random permutation, and use it to permute stratify
        final int[] permutation = getRange(stratify.length);
        final int[] permuted_stratify = slice(stratify, permutation);
        2. find the indices that would sort the permuted stratify array. The permutation ensures that entries with
           the same value in stratify will be in random order. However, they point to the wrong stratify values.
        final int[] permuted_sort_indices = argsort(permuted_stratify);
        3. unpermute the sort_indices to permute to the correct stratify values. Because the sort visited the values
           in random order, the ordering of stratify indices will be permuted (between indices pointing to equal
           stratify values).
        final int[] stratify_inds = slice(permutation, permuted_sort_indices);
        */
        final int[] permutation = MachineLearningUtils.getRandomPermutation(random, stratify.length);
        return MachineLearningUtils.slice(permutation, MachineLearningUtils.argsort(MachineLearningUtils.slice(stratify, permutation)));
    }

    public int[] getTrainRows() {
        return trainRows;
    }

    public int[] getTestRows() {
        return testRows;
    }

    private static class FoldSplitIterator implements Iterator<TrainTestSplit> {
        private final int[] split_index_ordering;
        private final int numFolds;
        private int fold;

        FoldSplitIterator(final int[] split_index_ordering, final int numFolds) {
            this.split_index_ordering = split_index_ordering;
            this.numFolds = numFolds;
            this.fold = 0;
        }

        @Override
        public boolean hasNext() {
            return fold < numFolds;
        }

        @Override
        public TrainTestSplit next() {
            final int numTest = 1 + (split_index_ordering.length - 1 - fold) / numFolds;
            final int numTrain = split_index_ordering.length - numTest;
            int[] testRows = new int[numTest];
            int[] trainRows = new int[numTrain];
            int trainIndex;
            int orderingIndex;
            if(fold > 0) {
                for(trainIndex = 0; trainIndex < fold; ++trainIndex) {
                    trainRows[trainIndex] = split_index_ordering[trainIndex];
                }
                orderingIndex = trainIndex;
            } else {
                orderingIndex = 0;
                trainIndex = 0;
            }
            for(int testIndex = 0; testIndex < testRows.length; ++testIndex) {
                testRows[testIndex] = split_index_ordering[orderingIndex];
                final int orderingStop = Math.min(orderingIndex + numFolds, split_index_ordering.length);
                for(++orderingIndex; orderingIndex < orderingStop; ++orderingIndex, ++trainIndex) {
                    trainRows[trainIndex] = split_index_ordering[orderingIndex];
                }
            }

            ++fold;
            Arrays.sort(trainRows);
            Arrays.sort(testRows);
            return new TrainTestSplit(trainRows, testRows);
        }
    }
}
