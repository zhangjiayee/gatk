package org.broadinstitute.hellbender.tools.spark.sv.ml;

import org.apache.logging.log4j.Logger;

import java.util.Random;

import static org.broadinstitute.hellbender.tools.spark.sv.ml.ClassifierTuner.DEFAULT_NUM_STRATIFY_BINS;

// we most likely need this split anyway, so I think it's useful to create one, not here, but in GATKXGBoostConstants or its own 1st level class
public class TrainValidateTestConfiguration {
    private final TruthSet tuneSet;
    private final int[] tuneStratify;
    private final TruthSet validateSet;
    private final int[] validateStratify;

    public TrainValidateTestConfiguration(final TruthSet truthSet,
                                          final Double hyperparameterTuningProportion,
                                          final Random random,
                                          final int numCrossvalidationFolds,
                                          final Logger localLogger) {

        // Choose min counts per stratify value so that after train-test split, each fold of cross-validation has at least
        // one element of each stratify value.
        final int minCountsPerStratifyValue = (hyperparameterTuningProportion == 0 || hyperparameterTuningProportion == 1) ?
                numCrossvalidationFolds
                : (int)Math.ceil(numCrossvalidationFolds / Math.min(hyperparameterTuningProportion, 1.0 - hyperparameterTuningProportion));
        localLogger.info("Stratifying data matrix to balance data splitting / cross-validation folds.");
        final int[] stratify = truthSet.getStratifyArray(
                DEFAULT_NUM_STRATIFY_BINS, minCountsPerStratifyValue
        );

        localLogger.info("Splitting data");
        final TrainTestSplit hyperSplit = TrainTestSplit.getTrainTestSplit(
                hyperparameterTuningProportion, truthSet.getNumRows(), random, stratify
        );
        tuneSet = truthSet.sliceRows(hyperSplit.getTrainRows());
        tuneStratify = MachineLearningUtils.slice(stratify, hyperSplit.getTrainRows());
        validateSet = truthSet.sliceRows(hyperSplit.getTestRows());
        validateStratify = MachineLearningUtils.slice(stratify, hyperSplit.getTestRows());
    }

    public TruthSet getTuneSet() {
        return tuneSet;
    }

    public int[] getTuneStratify() {
        return tuneStratify;
    }

    public TruthSet getValidateSet() {
        return validateSet;
    }

    public int[] getValidateStratify() {
        return validateStratify;
    }
}
