package org.broadinstitute.hellbender.tools.spark.sv.ml;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class GATKXGBoosterUnitTest extends GATKBaseTest {
    private static final String SV_UTILS_TEST_DIR = toolsTestDir + "spark/sv/utils/";
    private static final String TEST_MATRIX_DATA_FILE = SV_UTILS_TEST_DIR + "agaricus-integers.csv.gz";
    private static final TruthSet TEST_SET = TruthSet.loadCsvFile(TEST_MATRIX_DATA_FILE);
    private static final GATKXGBooster classifier = new GATKXGBooster();
    private static final int NUM_TRAINING_ROUNDS = GATKXGBoostConstants.DEFAULT_NUM_TRAINING_ROUNDS;
    private static final int NUM_EVAL_METRIC_TRAINING_ROUNDS = 100; // for tests where accurate results are not needed
    private static final int EARLY_STOPPING_ROUNDS = GATKXGBoostConstants.DEFAULT_EARLY_STOPPING_ROUNDS;
    private static final int NUM_CROSSVALIDATION_FOLDS = ClassifierTuner.DEFAULT_NUM_CROSSVALIDATION_FOLDS;
    private static final int NUM_TUNING_ROUNDS = ClassifierTuner.DEFAULT_NUM_TUNING_ROUNDS; // keep tests quick
    private static final double TUNING_FRACTION = 1.0 / (1.0 + NUM_TUNING_ROUNDS);
    private static final Map<String, Object> CLASSIFIER_PARAMS = GATKXGBoostConstants.DEFAULT_CLASSIFIER_PARAMETERS;
    private static final String[] TEST_EVAL_METRICS = {
            "rmse", "mae", "logloss", "error", "error@0.75", "auc", "ndcg", "map", "map@7000", "poisson-nloglik"
    };
    private static final Random random = new Random(GATKXGBoostConstants.DEFAULT_SEED);

    private static final double MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY = 0.99;

    private static void assertLabelsEqual(final int[] actuals, final int[] expecteds, final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], "at index=" + index + ": " + message);
        }
    }

    @Test(groups = "sv")
    protected void testBasicTrain() throws IOException {
        // check that a classifier can be trained from data
        final Map<String, Object> classifierParams = new HashMap<>(CLASSIFIER_PARAMS);
        classifierParams.put(ClassifierTuner.NUM_TRAINING_ROUNDS_KEY, NUM_EVAL_METRIC_TRAINING_ROUNDS);
        classifier.train(classifierParams, TEST_SET);

        // check that the classifier can make predictions on the same data;
        final int[] predictedLabels = classifier.predictClassLabels(TEST_SET.getFeatures());

        // check that the predictions match the data (this data is easy to predict, they should)
        assertLabelsEqual(predictedLabels, TEST_SET.getClassLabels(),
                "Predicted labels not identical to actual labels");

        // predict probabilities of whole matrix
        final double[][] probabilities = classifier.predictProbability(TEST_SET.getFeatures());
        // check that you get indentical results when predicting row-by-row
        for(int row = 0; row < TEST_SET.getNumRows(); ++row) {
            final double[] rowProbability = classifier.predictProbability(TEST_SET.getFeatures().getRow(row));
            assertArrayEquals(rowProbability, probabilities[row], 0,
                    "Row " + row + ": different probabilities predicted for matrix and row-by-row");
        }

        // save classifier to temporary file
        File tempFile = File.createTempFile("gatk-xgboost-classifier", "kryo");
        classifier.save(tempFile.getAbsolutePath());
        // load classifier from temporary file

        final GATKClassifier loadedClassifier = GATKClassifier.load(
                tempFile.getAbsolutePath()
        );
        // check that you get identical results to first probability predictions
        final double[][] loadedProbabilities = loadedClassifier.predictProbability(TEST_SET.getFeatures());
        assertMatrixEquals(loadedProbabilities, probabilities, 0.0,
                "Probabilities predicted by loaded classifier not equal to original");
    }

    @Test(groups = "sv")
    protected void testGetTrainingTrace() {
        final int[] stratify = TEST_SET.getClassLabels();

        final TrainTestSplit hyperSplit = TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_SET.getNumRows(), random, stratify
        );
        final TruthSet trainSet = TEST_SET.sliceRows(hyperSplit.getTrainRows());
        final TruthSet validateSet = TEST_SET.sliceRows(hyperSplit.getTestRows());

        final double[] trainingTrace = classifier.trainAndReturnQualityTrace(CLASSIFIER_PARAMS, trainSet, validateSet,
                NUM_TRAINING_ROUNDS, EARLY_STOPPING_ROUNDS);
        Assert.assertEquals(
                trainingTrace.length, NUM_TRAINING_ROUNDS,
                "Training trace did not have requested number of rounds (" + trainingTrace.length
                        + " instead of " + NUM_TRAINING_ROUNDS
        );
    }

    @Test(groups = "sv")
    protected void testGetMaximizeEvalMetric() {
        final int[] stratify = TEST_SET.getClassLabels();

        // Just get a small subset of data, and measure what the training algorithm is trying to do on that data.
        final TrainTestSplit hyperSplit = TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_SET.getNumRows(), random, stratify
        );
        final TruthSet trainSet = TEST_SET.sliceRows(hyperSplit.getTrainRows());
        final TruthSet validateSet = trainSet;

        final GATKXGBooster classifier = new GATKXGBooster();
        final Map<String, Object> classifierParams = new HashMap<>(GATKXGBoostConstants.DEFAULT_CLASSIFIER_PARAMETERS);
        for(final String evalMatric : TEST_EVAL_METRICS) {
            classifierParams.put(GATKXGBoostConstants.EVAL_METRIC_KEY, evalMatric);
            final boolean maximizeEvalMetric = classifier.getMaximizeEvalMetric(classifierParams);
            final double[] trainingTrace = classifier.trainAndReturnQualityTrace(
                    classifierParams, trainSet, validateSet, NUM_EVAL_METRIC_TRAINING_ROUNDS, Integer.MAX_VALUE);
            final double startVal = trainingTrace[0];
            final DoubleSummaryStatistics traceStats = Arrays.stream(trainingTrace).summaryStatistics();
            final double maxVal = traceStats.getMax();
            final double minVal = traceStats.getMin();
            final boolean traceIncreasedMoreThanDecreased = maxVal - startVal > startVal - minVal;
            Assert.assertEquals(traceIncreasedMoreThanDecreased, maximizeEvalMetric,
                    "For " + GATKXGBoostConstants.EVAL_METRIC_KEY + "=" + evalMatric + ": maximizeEvalMetric="
                            + maximizeEvalMetric + ", but traceIncreasedMoreThanDecreased=" + traceIncreasedMoreThanDecreased);
        }
    }

    @Test(groups = "sv")
    protected void testCrossvalidatedTuneAndTrain() {
        final int minCountsPerStratifyValue = (int)Math.round(NUM_CROSSVALIDATION_FOLDS / TUNING_FRACTION);
        final int[] stratify = TEST_SET.getStratifyArray(ClassifierTuner.DEFAULT_NUM_STRATIFY_BINS,
                                                         minCountsPerStratifyValue);

        final TrainTestSplit hyperSplit = TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_SET.getNumRows(), random, stratify
        );
        final TruthSet tuneSet = TEST_SET.sliceRows(hyperSplit.getTrainRows());
        final int[] tuneStratify = MachineLearningUtils.slice(stratify, hyperSplit.getTrainRows());
        final TruthSet validateSet =TEST_SET.sliceRows(hyperSplit.getTestRows());
        final int[] validateStratify = MachineLearningUtils.slice(stratify, hyperSplit.getTestRows());

        // set the number of threads
        final int chosenNumThreads = classifier.chooseNumThreads(CLASSIFIER_PARAMS, tuneSet, GATKXGBoostConstants.NUM_THREADS_KEY, logger);
        CLASSIFIER_PARAMS.put(GATKXGBoostConstants.NUM_THREADS_KEY, chosenNumThreads);

        final Map<String, Object> bestClassifierParameters = ClassifierTuner.tuneClassifierParameters(
                classifier,
                CLASSIFIER_PARAMS, GATKXGBoostConstants.DEFAULT_TUNING_PARAMETERS,
                ClassifierTuner.ClassifierTuningStrategy.RANDOM,
                tuneSet, random, tuneStratify, NUM_CROSSVALIDATION_FOLDS, NUM_TRAINING_ROUNDS,
                EARLY_STOPPING_ROUNDS, NUM_TUNING_ROUNDS
        );

        final int[] predictedTestLabels = classifier.crossvalidatePredict(
                validateSet, bestClassifierParameters, random, validateStratify, NUM_CROSSVALIDATION_FOLDS
        );

        final double accuracy = MachineLearningUtils.getPredictionAccuracy(
                predictedTestLabels, validateSet.getClassLabels()
        );

        Assert.assertTrue(accuracy >= MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY,
                "Crossvalidated prediction accuracy (" + accuracy + ") less than passing ("
                        + MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY + ")");
    }

    private static void assertArrayEquals(final double[] actuals, final double[] expecteds, final double tol,
                                          final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], tol, "at index=" + index + ": " + message);
        }
    }

    private static void assertMatrixEquals(final double[][] actuals, final double[][] expecteds, final double tol,
                                           final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Number of rows not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            assertArrayEquals(actuals[index], expecteds[index], tol, "at row=" + index + ": " + message);
        }
    }
}
