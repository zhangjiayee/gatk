package org.broadinstitute.hellbender.tools.spark.sv.ml;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;


/**
 * <p>(Internal) Demo of GATKXGBoostConstants + MachineLearningUtils: trains classifier on data in .csv file and saves classifier to binary file</p>
 *
 * <p><h3>Inputs</h3>
 * <ul>
 *     <li>A csv file of numeric data, with first column being class labels.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A binary file that stores a trained classifier.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk ExampleTrainXGBoostClassifier \
 *     -I input_data.csv \
 *     -O classifier.model
 * </pre></p>
 *
 * <p><h3>Caveats</h3>
 * <li>This tool uses xgboost's multi-threading. The user can pass --nthread to set the number of threads.</li>
 * <li>If --auto-select-nthread is passed, then nthread sets the upper bound on threads to use, and the actual value
 * will be obtained by finding the best performance when repeatedly solving a small test problem using varying numbers
 * of threads.</li>
 * <li>Because this tool is multi-threaded, it must be used on a single computer (not a multi-computer cluster).
 * However, the classifiers are KryoSerializable, so trained classifiers can be used to predict probabilities or class
 * labels in a spark environment.</li>
 * Coverage much lower than that probably won't work well.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Demo of GATKXGBoostConstants + MachineLearningUtils: trains classifier on data in .csv file and saves classifier to binary file",
        summary = "This tool tunes hyperparameters, assesses classifier quality, and trains a classifier to predict class" +
                        " membership from data provided in a .csv file. The trained classifier is saved to a binary model file.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class ExampleTrainXGBoostClassifier extends CommandLineProgram {
    private static final long serialVersionUID = 1L;
    private static final Logger localLogger = LogManager.getLogger(ExampleTrainXGBoostClassifier.class);

    @Argument(doc = "path to .csv file storing data to train classifier. It is expected that the file will have all numeric" +
            "values, and the first column will contain class label.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    private String demoDataFile;

    @Argument(doc = "full path to save output (binary) classifier model file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    private String classifierModelFile;

    @Argument(doc="Stop classifier training if score does not improve for this many consecutive rounds.",
            fullName = "early-stopping-rounds", optional = true)
    private final int earlyStoppingRounds = GATKXGBoostConstants.DEFAULT_EARLY_STOPPING_ROUNDS;

    @Argument(doc="Train classifier for at most this many rounds.",
            fullName = "max-training-rounds", optional = true)
    private final int maxTrainingRounds = GATKXGBoostConstants.DEFAULT_NUM_TRAINING_ROUNDS;

    @Argument(doc="When performing cross-validation, use this many folds.",
            fullName = "num-crossvalidation-folds", optional = true)
    private final int numCrossvalidationFolds = ClassifierTuner.DEFAULT_NUM_CROSSVALIDATION_FOLDS;

    @Argument(doc="When optimizing hyperparameters, search for this many rounds.",
            fullName = "num-hyperparameter-optimization-rounds", optional = true)
    private final int numTuningRounds = ClassifierTuner.DEFAULT_NUM_TUNING_ROUNDS;

    @Argument(doc="When optimizing hyperparameters, reserve this proportion of data for tuning hyperparameters.",
            fullName = "hyperparameter-tuning-proportion", optional = true)
    private Double hyperparameterTuningProportion = null;

    @Argument(doc="Use this metric to evaluate performance of the classifier.",
            fullName = "eval-metric", optional = true)
    private final String evalMetric = GATKXGBoostConstants.DEFAULT_EVAL_METRIC;

    @Argument(doc="Seed for random numbers. If null, initialize randomly",
            fullName = "random-seed", optional = true)
    private final Long seed = GATKXGBoostConstants.DEFAULT_SEED;

    @Advanced
    @Argument(doc="Auto-select optimal number of threads to use for training classifier? If false, use nthread as"
                  + " specified by user. If false, interpret that number as an upper bound, and select number of threads"
                  + " with maximal throughput on small test problem",
              fullName = "auto-select-nthread", optional = true)
    private final Boolean autoSelectNumThreads = false;

    @Argument(doc="Number of threads to use for training classifier. If ",
            fullName = "nthread", optional = true)
    private final int numThreads = autoSelectNumThreads ? Runtime.getRuntime().availableProcessors()
            : (int)Math.round(GATKXGBoostConstants.GUESS_OPTIMAL_NUM_THREADS_PROPORTION * Runtime.getRuntime().availableProcessors());

    @Argument(doc="Tuning strategy for choosing classifier hyperparameters",
            fullName = "classifier-tuning-strategy", optional = true)
    final ClassifierTuner.ClassifierTuningStrategy classifierTuningStrategy
            = ClassifierTuner.ClassifierTuningStrategy.RANDOM;

    private final Random random = (seed == null ? new Random() : new Random(seed));
    /**
     * Demo stuff ENDS here....
     */

    @Override
    protected void onStartup() {
        if(hyperparameterTuningProportion == null) {
            // This proportion produces the fastest training times overall, not necessarily the most useful results.
            hyperparameterTuningProportion = 1.0 / (1.0 + numTuningRounds);
        }
    }

    @Override
    protected Object doWork() {
        localLogger.info("Loading demo data");
        final TruthSet truthSet = TruthSet.loadCsvFile(demoDataFile);

        // tune hyperparameters, assess cross-validated classifier accuracy, train final classifier
        final GATKClassifier classifier = trainClassifierAndAssessAccuracy(truthSet);

        localLogger.info("Saving final classifier to " + classifierModelFile);
        try {
            classifier.save(classifierModelFile);
        } catch(IOException err) {
            throw new GATKException(err.getClass() + ": " + err.getMessage());
        }

        // load trained classifier from disk, and use it to predict class membership probabilities, or class labels
        loadClassifierAndPredictClass(truthSet.getFeatures());

        return classifier;
    }

    /**
     * Demo tuning hyperparameters, assessing classifier accuracy, and training final classifier
     */
    private GATKClassifier trainClassifierAndAssessAccuracy(final TruthSet truthSet) {

        final TrainValidateTestConfiguration trainValidateTestConfiguration =
                new TrainValidateTestConfiguration(truthSet, hyperparameterTuningProportion, random, numCrossvalidationFolds, localLogger);

        final Map<String, Object> bestClassifierParameters = getBestClassifierParameters(trainValidateTestConfiguration);

        localLogger.info("Training final classifier");
        final GATKXGBooster classifier = new GATKXGBooster();
        classifier.train(bestClassifierParameters, truthSet);
        return classifier;
    }

    private Map<String, Object> getBestClassifierParameters(final TrainValidateTestConfiguration trainValidateTestConfiguration) {

        localLogger.info("Creating classifier object");
        final GATKXGBooster gatkxgBooster = new GATKXGBooster();

        final TruthSet tuneSet = trainValidateTestConfiguration.getTuneSet();
        final int[] tuneStratify = trainValidateTestConfiguration.getTuneStratify();
        final TruthSet validateSet = trainValidateTestConfiguration.getValidateSet();
        final int[] validateStratify = trainValidateTestConfiguration.getValidateStratify();

        Map<String, Object> bestClassifierParameters = new HashMap<>(GATKXGBoostConstants.DEFAULT_CLASSIFIER_PARAMETERS);
        bestClassifierParameters.put(GATKXGBoostConstants.EVAL_METRIC_KEY, evalMetric);

        // set the number of threads
        if(autoSelectNumThreads) {
            final int chosenNumThreads = gatkxgBooster.chooseNumThreads(
                    bestClassifierParameters, tuneSet, GATKXGBoostConstants.NUM_THREADS_KEY, localLogger);
            bestClassifierParameters.put(GATKXGBoostConstants.NUM_THREADS_KEY, chosenNumThreads);
            localLogger.info("Chose " + chosenNumThreads + " threads");
        } else {
            bestClassifierParameters.put(GATKXGBoostConstants.NUM_THREADS_KEY, numThreads);
            localLogger.info("User selected " + numThreads + " threads");
        }

        // Note: tuning hyperparameters is basically always necessary: if data has changed enough to merit training a
        // new classifier, then it merits finding the new best hyperparameters
        if(tuneSet.getNumRows() > 0) {
            localLogger.info("Tuning hyperparameters");
            bestClassifierParameters = ClassifierTuner.tuneClassifierParameters(
                    gatkxgBooster,
                    bestClassifierParameters, GATKXGBoostConstants.DEFAULT_TUNING_PARAMETERS, classifierTuningStrategy,
                    tuneSet, random, tuneStratify, numCrossvalidationFolds, maxTrainingRounds, earlyStoppingRounds,
                    numTuningRounds
            );
            localLogger.info("bestClassifierParameters: " + bestClassifierParameters.toString());
        } else {
            localLogger.info("skipping tuning hyperparameters, using default parameters");
        }

        // note cross-validation is necessary if you want to estimate the accuracy of the new classifier. If you don't
        // care (e.g. are committed to using the classifier regardless) then you can set hyperparameterTuningProportion
        // to 1.0, or skip splitting the data and just tune on the whole set.
        if(validateSet.getNumRows() > 0) {
            localLogger.info("Cross-val predicting");
            final int[] predictedTestLabels = gatkxgBooster.crossvalidatePredict(
                    validateSet, bestClassifierParameters, random, validateStratify, numCrossvalidationFolds
            );

            final double accuracy = MachineLearningUtils.getPredictionAccuracy(
                    predictedTestLabels, validateSet.getClassLabels()
            );
            localLogger.info("Crossvalidated accuracy = " + String.format("%.1f%%", 100.0 * accuracy));
        } else {
            localLogger.info("skipping evaluating crossvalidated accuracy");
        }

        return bestClassifierParameters;
    }

    /**
     * Demo loading saved classifier and using it to make predictions on data.
     */
    private void loadClassifierAndPredictClass(final RealMatrix dataMatrix) {
        localLogger.info("Re-loading saved classifier");
        final GATKClassifier loadedClassifier;
        try {
            loadedClassifier = GATKClassifier.load(classifierModelFile);
        } catch(IOException err) {
            throw new GATKException(err.getClass() +": " + err.getMessage());
        }

        // These aren't used here, but presumably any real application would do something based on probability or label
        localLogger.info("predicting probability that data is in each class");
        final double[][] probabilities = loadedClassifier.predictProbability(dataMatrix);

        localLogger.info("predicting class labels");
        final int[] classLabels = loadedClassifier.predictClassLabels(dataMatrix);
    }
}
