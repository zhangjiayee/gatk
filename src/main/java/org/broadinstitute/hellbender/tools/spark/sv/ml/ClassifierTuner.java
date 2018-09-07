package org.broadinstitute.hellbender.tools.spark.sv.ml;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Abstract base class for tuning classifier hyperparameters.
 * The base class manages keeping track of previously
 * seen / best hyperparameters and corresponding scores.
 * Actual classifier tuners must extend ClassifierTuner
 * and override {@link #chooseNextHyperparameters()}.
 * Note, they will probably always be called by
 * {@link #tuneClassifierParameters(GATKClassifier, Map, Map, ClassifierTuningStrategy, TruthSet, Random, int[], int, int, int, int)}
 * rather than instantiated directly.
 */
abstract class ClassifierTuner {
    public static final int DEFAULT_NUM_CROSSVALIDATION_FOLDS = 5;
    public static final int DEFAULT_NUM_TUNING_ROUNDS = 100;
    public static final int DEFAULT_NUM_STRATIFY_BINS = 5;
    public static final String NUM_TRAINING_ROUNDS_KEY = "num_training_rounds";

    // To-do: write ClassifierTuner with strategy = BAYES
    public enum ClassifierTuningStrategy { RANDOM }

    protected final GATKClassifier classifier;
    protected final Map<String, Object> classifierParameters;
    protected final Map<String, ClassifierParamRange<?>> tuneParameters;
    protected final List<Map<String, Object>> hyperparameterSets;
    protected final List<Double> hyperparameterScores;
    protected final boolean maximizeEvalMetric;

    // keep these as class members instead of locals so that chooseNextHyperparameters() can see them if it needs to
    protected Map<String, Object> bestParameters;
    protected double bestScore;
    protected int numTuningRounds;

    private final TruthSet truthSet;
    private final List<TrainTestSplit> splits;
    private final int maxTrainingRounds;
    private final int earlyStoppingRounds;

    ClassifierTuner(final GATKClassifier classifier, final Map<String, Object> classifierParameters,
                    final Map<String, ClassifierParamRange<?>> tuneParameters, final TruthSet truthSet,
                    final List<TrainTestSplit> splits, final int maxTrainingRounds, final int earlyStoppingRounds,
                    final int numTuningRounds) {
        if(numTuningRounds < 1) {
            throw new IllegalArgumentException("numTuningRounds (" + numTuningRounds + ") must be >= 1");
        }

        this.classifier = classifier;
        this.tuneParameters = tuneParameters;
        this.hyperparameterSets = new ArrayList<>();
        this.hyperparameterScores = new ArrayList<>();
        this.classifierParameters = classifierParameters;

        this.truthSet = truthSet;
        this.splits = splits;
        this.maxTrainingRounds = maxTrainingRounds;
        this.earlyStoppingRounds = earlyStoppingRounds;
        this.numTuningRounds = numTuningRounds;

        maximizeEvalMetric = classifier.getMaximizeEvalMetric(classifierParameters);
        bestParameters = null;
        bestScore = maximizeEvalMetric ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
    }

    /**
     * Return next candidate set of hyperparameters
     */
    abstract protected Map<String, Object> chooseNextHyperparameters();

    // TODO: 9/7/18 is this for tuning hyper-parameters or model parameters?
    /**
     * Obtain optimal classifierParameters.
     * Get quality estimates by calling trainAndReturnQualityTrace() with
     * cross-validation.
     * Select optimal number of training rounds (to be used when early stopping is not possible).
     * @param classifier           classifier whose parameters to be tuned??????
     * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
     *                             initialization and/or training routine. Since map values are Objects, the
     *                             classifier must keep track of their type. This must specify any parameters that
     *                             will not be tuned. If there is overlap in keys between classifierParameters and
     *                             tuneClassifierParameters, the tuned value will take precedence.
     * @param tuneClassifierParameters Map (from parameter names to ClassifierParamRange) specifying values that
     *                                 will be tuned, their range of allowed values, and how they are distributed
     *                                 (float vs integer, linear vs log).
     * @param classifierTuningStrategy enum specifying strategy for selecting candidate parameter values
     * @param truthSet                 Data with matrix of features paired to correct class labels.
     * @param random                   Random number generator.
     * @param stratify                 Vector of integers listing stratification class. Cross-validation will split
     *                                 data so that each stratify value has balanced number of instances in each
     *                                 fold. If null, split randomly. It is highly desirable to stratify at least
     *                                 based on class label.
     * @param numCrossvalidationFolds  Must be >= 2
     * @param maxTrainingRounds        Train classifiers with candidate hyperparameters for at most this many rounds.
     * @param earlyStoppingRounds      Stop training if the best training quality was this many rounds ago.
     * @param numTuningRounds          Try this many candidate hyperparameters before stopping and selecting the best.
     * @return bestHyperparameters     classifierParameters updated with optimal values from tuneClassifierParameters.
     */
    static Map<String, Object> tuneClassifierParameters(final GATKClassifier classifier,
                                                        final Map<String, Object> classifierParameters,
                                                        final Map<String, ClassifierParamRange<?>> tuneClassifierParameters,
                                                        final ClassifierTuningStrategy classifierTuningStrategy,
                                                        final TruthSet truthSet, final Random random,
                                                        final int[] stratify, final int numCrossvalidationFolds,
                                                        final int maxTrainingRounds, final int earlyStoppingRounds,
                                                        final int numTuningRounds) {

        final ClassifierTuner classifierTuner = getClassifierTuner(classifier, truthSet, random,
                classifierParameters, tuneClassifierParameters, classifierTuningStrategy,
                stratify, numCrossvalidationFolds, maxTrainingRounds, earlyStoppingRounds, numTuningRounds);

        return classifierTuner.getBestHyperParameters();
    }

    // TODO: 9/7/18 not sure if I understand the intention of the original code
    /**
     * Given various input parameters, construct a {@link ClassifierTuner}.
     */
    private static ClassifierTuner getClassifierTuner(final GATKClassifier classifier,
                                                      final TruthSet truthSet,
                                                      final Random random,
                                                      final Map<String, Object> classifierParameters,
                                                      final Map<String, ClassifierParamRange<?>> tuneClassifierParameters,
                                                      final ClassifierTuningStrategy classifierTuningStrategy,
                                                      final int[] stratify,
                                                      final int numCrossvalidationFolds,
                                                      final int maxTrainingRounds,
                                                      final int earlyStoppingRounds,
                                                      final int numTuningRounds) {

        final List<TrainTestSplit> splits = new ArrayList<>(numCrossvalidationFolds);
        TrainTestSplit.getCrossvalidationSplits(
                numCrossvalidationFolds, truthSet.getNumRows(), random, stratify
        ).forEachRemaining(splits::add);

        final ClassifierTuner classifierTuner;
        switch(classifierTuningStrategy) {
            case RANDOM:
                classifierTuner = new RandomClassifierTuner(
                        classifier, classifierParameters, tuneClassifierParameters, truthSet,
                        splits, maxTrainingRounds, earlyStoppingRounds, numTuningRounds, random
                );
                break;
            default:
                throw new IllegalStateException("Invalid ClassifierTuningStrategy: " + classifierTuningStrategy);
        }
        return classifierTuner;
    }

    /**
     * Try multiple candidate hyperparameters (chosen according to strategy in non-abstract ClassifierTuner derived
     * class). Keep track of best, and after numTuningRounds candidates have been tried, return the best hyperparameters
     */
    private Map<String, Object> getBestHyperParameters() {
//        MachineLearningUtils.localLogger.info("Getting best parameters"); // I'm removing this because I think the caller of this method should be responsible for logging steps.
        ClassifierParamRange.ConsoleProgressBar progress = new ClassifierParamRange.ConsoleProgressBar(numTuningRounds);
        for(int i = 0; i < numTuningRounds; ++i) {
            final Map<String, Object> hyperparameters = chooseNextHyperparameters();
            hyperparameterSets.add(hyperparameters);
            final Map<String, Object> testParameters = new HashMap<>(classifierParameters);
            testParameters.putAll(hyperparameters);
            final double[][] trainingTraces = classifier.getCrossvalidatedTrainingTraces(
                    testParameters, truthSet, splits, maxTrainingRounds, earlyStoppingRounds
            );
            final double score = classifier.getTrainingScore(testParameters, trainingTraces);
            hyperparameterScores.add(score);
            // This is the new best score if a) it is better than the previous best score so far -OR-
            //                               b) it exactly ties the best score, but uses fewer rounds of training
            if((maximizeEvalMetric ? score > bestScore : score < bestScore)
                    || (score == bestScore &&
                    (int)testParameters.get(NUM_TRAINING_ROUNDS_KEY) < (int)bestParameters.get(NUM_TRAINING_ROUNDS_KEY))) {
                bestScore = score;
                bestParameters = testParameters;
            }
            progress.update(1);
        }
        return bestParameters;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Concrete ClassifierTuner that selects hyperparameters randomly.
     * Specifically:
     * 1) each parameter value is sampled evenly (log or linear) along its range, with numTuningRounds samples.
     * 2) these samples are permuted into random order
     * 3) different parameter values are permuted independently
     * The result is an even sampling in the hyper-rectangle, but with less clumping.
     * Unimportant parameters (e.g. those that don't affect quality) will not decrease
     * sample density in important parameters.
     * Thus this method should outperform grid search with an equivalent {@link #numTuningRounds}.
     */
    static class RandomClassifierTuner extends ClassifierTuner {
        private final Map<String, Object[]> randomParameters;

        RandomClassifierTuner(final GATKClassifier classifier,
                              final Map<String, Object> classifierParameters,
                              final Map<String, ClassifierParamRange<?>> tuneParameters,
                              final TruthSet truthSet,
                              final List<TrainTestSplit> splits,
                              final int maxTrainingRounds,
                              final int earlyStoppingRounds,
                              final int numTuningRounds,
                              final Random random) {
            super(classifier, classifierParameters, tuneParameters, truthSet, splits, maxTrainingRounds,
                    earlyStoppingRounds, numTuningRounds);
            randomParameters = tuneParameters.entrySet().stream().collect(
                    Collectors.toMap(
                            Map.Entry::getKey,
                            x -> x.getValue().getRandomSamples(random, numTuningRounds)
                    )
            );
        }

        @Override
        protected Map<String, Object> chooseNextHyperparameters() {
            final int index = hyperparameterScores.size();
            return randomParameters.entrySet().stream().collect(
                    Collectors.toMap(
                            Map.Entry::getKey,
                            x -> x.getValue()[index]
                    )
            );
        }
    }
}
