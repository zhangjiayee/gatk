package org.broadinstitute.hellbender.tools.spark.sv.ml;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.KryoSerializable;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.Logger;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;

/**
 * Abstract base classifier class. Actual classifiers must extend GATKClassifier and wrap whatever 3rd-party package
 * implements the actual numerics. Note that implementing class must be Serializable and KryoSerializble (that is
 * how saving is implemented)
 */
public abstract class GATKClassifier implements Serializable, KryoSerializable {
    private static final long serialVersionUID = 1L;
    private static final int NUM_CALIBRATION_TRAINING_ROUNDS = 500;
    private static final int NUM_CALIBRATION_CLASS_ROWS = 500;

    // TODO: 9/7/18 I guess the reason why parameters are not part of this object is due to serialization?
    private final double[][] singlePredictWrapper = new double [1][];

    /* public abstract routines. Implementing classes must override these */

    /**
     * Train classifier.
     * This routine must be overridden by implementing class. It should update "this" to a trained state
     * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
     *                             initialization and/or training routine. Since map values are Objects, the
     *                             classifier must keep track of their type.
     * @param truthSet             Data with matrix of features paired to correct class labels.
     * @return upon success, the function should return "this"
     */
    public abstract GATKClassifier train(final Map<String, Object> classifierParameters,
                                         final TruthSet truthSet);

    /**
     * Train classifier and return trace of quality vs training round.
     * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
     *                             initialization and/or training routine. Since map values are Objects, the
     *                             classifier must keep track of their type.
     * @param trainingSet          Truth data used for updating classifier during training.
     * @param evaluationSet        Truth data used for evaluating classifier performance during training.
     * @param maxTrainingRounds    Train for no more than this many rounds.
     * @param earlyStoppingRounds  Training should stop if the best quality is more than earlyStoppingRounds ago.
     *                             If early stopping is triggered, remaining trace values should be set to last
     *                             obtained training value (i.e. not the best).
     * @return double[] with length maxTrainingRounds, each value storing classifier "quality" as a function of
     *         training round.
     */
    public abstract double[] trainAndReturnQualityTrace(
            final Map<String, Object> classifierParameters, final TruthSet trainingSet,
            final TruthSet evaluationSet, final int maxTrainingRounds, final int earlyStoppingRounds);

    /**
     * return numDataRows x numClasses double[][] with probability each data point is a member of each class
     */
    public abstract double[][] predictProbability(final RealMatrix matrix);

    /**
     * return true if the "quality" returned by trainAndReturnQualityTrace should be maximized, false if it should
     * be minimized
     */
    public abstract boolean getMaximizeEvalMetric(final Map<String, Object> classifierParameters);

    /* public defined methods. Can be overridden if more efficient routines are available for specific implementation */

    /**
     * Predict the probability that an individual data point is a member of each class
     */
    public double[] predictProbability(final double[] featureVector) {
        singlePredictWrapper[0] = featureVector;
        final RealMatrix matrix = new Array2DRowRealMatrix(singlePredictWrapper, false);
        return (predictProbability(matrix)[0]);
    }

    /**
     * Predict the class label for each data point (represented by matrix rows)
     */
    public int[] predictClassLabels(final RealMatrix matrix) {
        final double [][] predictedProbabilities = predictProbability(matrix);
        final int [] predictedLabels = new int[predictedProbabilities.length];
        if(predictedProbabilities.length == 0) {
            return predictedLabels;
        }
        final int numColumns = predictedProbabilities[0].length;
        if(numColumns == 1) {
            // binary classifier, reporting only probability of class == 1 (or "true")
            for(int row = 0; row < predictedProbabilities.length; ++row) {
                predictedLabels[row] = predictedProbabilities[row][0] >= 0.5 ? 1 : 0;
            }

        } else {
            // multiclass classifier (or at binary independently reporting probability of class 0 or 1)
            for (int row = 0; row < predictedProbabilities.length; ++row) {
                predictedLabels[row] = MachineLearningUtils.argmax(predictedProbabilities[row]);
            }
        }
        return predictedLabels;
    }

    /**
     * Save classifier to file at specified path (via Kryo)
     */
    public void save(final String saveFilePath) throws IOException {
        try(FileOutputStream fileOutputStream = new FileOutputStream(saveFilePath)) {
            save(fileOutputStream);
        }
    }

    /**
     * Save classifier to specified FileOutputStream (via Kryo)
     */
    public void save(final FileOutputStream fileOutputStream) {
        final Kryo kryo = new Kryo();
        final Output output = new Output(fileOutputStream);
        kryo.writeClassAndObject(output, this);
        output.close();
    }

    /**
     * Load classifier from file at specified path (via Kryo)
     */
    public static GATKClassifier load(final String saveFilePath) throws IOException {
        try(FileInputStream fileInputStream = new FileInputStream(saveFilePath)) {
            return load(fileInputStream);
        }
    }

    /**
     * Load classifier from specified FileInputStream (via Kryo)
     */
    public static GATKClassifier load(final FileInputStream fileInputStream) {
        final Kryo kryo = new Kryo();
        final Input input = new Input(fileInputStream);
        final GATKClassifier classifier = (GATKClassifier)kryo.readClassAndObject(input);
        input.close();
        return classifier;
    }

    /**
     * Predict class labels using cross-validation, so classifier does not predict on data that it was trained on.
     * @param truthSet             Data with matrix of features paired to correct class labels.
     * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
     *                             initialization and/or training routine. Since map values are Objects, the
     *                             classifier must keep track of their type. This must specify any parameters that
     *                             will not be tuned. If there is overlap in keys between classifierParameters and
     *                             tuneClassifierParameters, the tuned value will take precedence.
     * @param random               Random number generator
     * @param stratify             Vector of integers listing stratification class. Cross-validation will split data
     *                             so that each stratify value has balanced number of instances in each fold. If
     *                             null, split randomly. It is highly desirable to stratify at least based on class
     *                             label.
     * @param numCrossvalidationFolds  >= 2
     */
    public int[] crossvalidatePredict(final TruthSet truthSet, final Map<String, Object> classifierParameters,
                                      final Random random, final int[] stratify, final int numCrossvalidationFolds) {
        final int[] predictedLabels = new int[truthSet.getNumRows()];
        final Iterator<TrainTestSplit> splitIterator = TrainTestSplit.getCrossvalidationSplits(
                numCrossvalidationFolds, truthSet.getNumRows(), random, stratify
        );
        while(splitIterator.hasNext()) {
            final TrainTestSplit split = splitIterator.next();
            // train on training data from this crossvalidation split
            train(classifierParameters, truthSet.sliceRows(split.getTrainRows()));
            // predict on testing data from this split
            final int[] predictedTestLabels = predictClassLabels(MachineLearningUtils.sliceRows(truthSet.getFeatures(), split.getTestRows()));
            // and assign those values into the final predictions
            MachineLearningUtils.sliceAssign(predictedLabels, split.getTestRows(), predictedTestLabels);
        }
        return predictedLabels;
    }

    /**
     * Choose the number of threads to obtain fastest training times.
     * Update classifierParameters with the appropriate value for the number of threads.
     * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
     *                             initialization and/or training routine. Since map values are Objects, the
     *                             classifier must keep track of their type. This must specify any parameters that
     *                             will not be tuned. If there is overlap in keys between classifierParameters and
     *                             tuneClassifierParameters, the tuned value will take precedence.
     * @param truthSet       Data to obtain training times from. Ideally this data should have similar structure to
     *                       actual data (there is no reason you can't subsequently train on it). Note that a small
     *                       subset of the data will be selected, so that training times are short and this routine
     *                       does not waste a lot of time.
     * @param numThreadsKey  Key name in classifierParameters for number of threads.
     * @return               chosen number of threads
     */
    public int chooseNumThreads(final Map<String, Object> classifierParameters, final TruthSet truthSet,
                                final String numThreadsKey, final Logger logger) {
        final int numCalibrationRows = NUM_CALIBRATION_CLASS_ROWS * 2;
        logger.info("numCalibrationRows = " + numCalibrationRows);
        logger.info("numTrainingRows = " + truthSet.getNumRows());

        final TruthSet calibrationMatrix = getCalibrationMatrix(truthSet, numCalibrationRows, logger);

        final Map<String, Object> calibrationParams = new HashMap<>(classifierParameters);
        calibrationParams.put(ClassifierTuner.NUM_TRAINING_ROUNDS_KEY, NUM_CALIBRATION_TRAINING_ROUNDS);
        final int maxNumThreads =
                classifierParameters.containsKey(numThreadsKey) && (int)classifierParameters.get(numThreadsKey) > 0 ?
                        (int)classifierParameters.get(numThreadsKey) : Runtime.getRuntime().availableProcessors();
        if(maxNumThreads == 1) {
//            classifierParameters.put(numThreadsKey, 1); personally, I try not to modify input arguments, but that's personal preference
            return 1;
        }

        return chooseNumbThreadsBasedOnTrainTime(numThreadsKey, calibrationMatrix, calibrationParams, maxNumThreads);
    }

    private int chooseNumbThreadsBasedOnTrainTime(final String numThreadsKey,
                                                  final TruthSet calibrationMatrix,
                                                  final Map<String, Object> calibrationParams,
                                                  final int maxNumThreads) {
        int chosenNumThreads = 1;
        long bestElapsedTime = Long.MAX_VALUE;
        for(int numThreads = 1; numThreads < maxNumThreads; ++numThreads) {
            calibrationParams.put(numThreadsKey, numThreads);
            final long elapsedTime = getTrainingTime(calibrationParams, calibrationMatrix);
            if(elapsedTime < bestElapsedTime) {
                bestElapsedTime = elapsedTime;
                chosenNumThreads = numThreads;
//                classifierParameters.put(numThreadsKey, numThreads); again personal preference
            }
        }
        return chosenNumThreads;
    }

    private static TruthSet getCalibrationMatrix(final TruthSet truthSet, final int numCalibrationRows, final Logger logger) {
        final TruthSet calibrationMatrix;
        if(truthSet.getNumRows() <= numCalibrationRows) {
            calibrationMatrix = truthSet;
        } else {
            final int[] stratify = truthSet.getClassLabels();
           logger.info("trainingFraction = " + numCalibrationRows / (double)truthSet.getNumRows());
            final TrainTestSplit trainTestSplit = TrainTestSplit.getTrainTestSplit(
                    numCalibrationRows / (double)truthSet.getNumRows(),
                    truthSet.getNumRows(), new Random(), stratify
            );
           logger.info("numTrain = " + trainTestSplit.getTrainRows().length);
           logger.info("numTest = " + trainTestSplit.getTestRows().length);
            calibrationMatrix = truthSet.sliceRows(trainTestSplit.getTrainRows());
        }
        return calibrationMatrix;
    }

    /* private methods */

    private long getTrainingTime(final Map<String, Object> classifierParameters, final TruthSet truthSet) {
        final long startTime = System.nanoTime();
        train(classifierParameters, truthSet);
        return System.nanoTime() - startTime;
    }

    double[][] getCrossvalidatedTrainingTraces(final Map<String, Object> classifierParameters,
                                               final TruthSet truthSet, final List<TrainTestSplit> splits,
                                               final int maxTrainingRounds, final int earlyStoppingRounds) {
        final int numCrossvalidationFolds = splits.size();
        double[][] trainingTraces = new double[numCrossvalidationFolds][];
        for(int fold = 0; fold < numCrossvalidationFolds; ++fold) {
            final TrainTestSplit split = splits.get(fold);
            trainingTraces[fold] = trainAndReturnQualityTrace(
                    classifierParameters, truthSet.sliceRows(split.getTrainRows()), truthSet.sliceRows(split.getTestRows()),
                    maxTrainingRounds, earlyStoppingRounds
            );
        }
        return trainingTraces;
    }

    double getTrainingScore(final Map<String, Object> classifierParameters,
                            final double[][] trainingTraces) {
        // find the index with the best total (i.e. mean) score across rounds. This yields the optimal number of
        // rounds of training
        final boolean maximizeEvalMetric = getMaximizeEvalMetric(classifierParameters);
        double bestTotalScore = maximizeEvalMetric ? Double.MIN_VALUE : Double.MAX_VALUE;
        int bestRoundIndex = -1;
        final int maxTrainingRounds = trainingTraces[0].length;
        for(int roundIndex = 0; roundIndex < maxTrainingRounds; ++roundIndex) {
            double roundScore = trainingTraces[0][roundIndex];
            for(int fold = 1; fold < trainingTraces.length; ++fold) {
                roundScore += trainingTraces[fold][roundIndex];
            }
            if(maximizeEvalMetric ? roundScore > bestTotalScore : roundScore < bestTotalScore) {
                // the score at this round of training is the best so far
                bestTotalScore = roundScore;
                bestRoundIndex = roundIndex;
            }
        }
        final int numTrainingRounds = bestRoundIndex + 1;

        // report the overall score for this set of parameters as the score of the *worst* trace at the optimal training
        // index. Selecting the worst trace demands high reliability from the classifier across similar data sets.
        double trainingScore = trainingTraces[0][bestRoundIndex];
        if(maximizeEvalMetric) {
            for (int fold = 1; fold < trainingTraces.length; ++fold) {
                trainingScore = Math.min(trainingScore, trainingTraces[fold][bestRoundIndex]);
            }
        } else {
            for (int fold = 1; fold < trainingTraces.length; ++fold) {
                trainingScore = Math.max(trainingScore, trainingTraces[fold][bestRoundIndex]);
            }
        }

        classifierParameters.put(ClassifierTuner.NUM_TRAINING_ROUNDS_KEY, numTrainingRounds);
        return trainingScore;
    }
}
