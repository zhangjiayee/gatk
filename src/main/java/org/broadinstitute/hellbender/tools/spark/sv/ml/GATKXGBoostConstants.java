package org.broadinstitute.hellbender.tools.spark.sv.ml;

import java.util.HashMap;
import java.util.Map;

/**
 * Concrete implementation of MachineLearningUtils that wraps XGBoost4J
 */
public class GATKXGBoostConstants extends MachineLearningUtils {
    public static final double GUESS_OPTIMAL_NUM_THREADS_PROPORTION = 0.5;

    public static final String LEARNING_RATE_KEY = "eta";
    public static final double DEFAULT_LEARNING_RATE = 0.5;
    public static final String MAX_DEPTH_KEY = "max_depth";
    public static final String GAMMA_KEY = "gamma";
    public static final String MIN_CHILD_WEIGHT_KEY = "min_child_weight";
    public static final String SUBSAMPLE_KEY = "subsample";
    public static final String COLSAMPLE_BY_TREE_KEY = "colsample_by_tree";
    public static final String COLSAMPLE_BY_LEVEL_KEY = "colsample_by_level";
    public static final String MAX_DELTA_STEP_KEY = "max_delta_step";

    public static final int DEFAULT_MAX_DEPTH = 6;
    public static final String SILENT_KEY = "silent";
    public static final int DEFAULT_SILENT = 1;
    public static final String OBJECTIVE_KEY = "objective";
    public static final String DEFAULT_OBJECTIVE_VALUE = "binary:logistic";
    public static final String EVAL_METRIC_KEY = "eval_metric";
    public static final String DEFAULT_EVAL_METRIC = "auc";
    public static final String NUM_THREADS_KEY = "nthread";
    public static final String SCALE_POS_WEIGHT_KEY = "scale_pos_weight";
    public static final String SEED_KEY = "seed";
    public static final long DEFAULT_SEED = 0L;
    public static final int DEFAULT_NUM_TRAINING_ROUNDS = 1000;
    public static final int DEFAULT_EARLY_STOPPING_ROUNDS = 50;
    @SuppressWarnings("serial")
    public static final Map<String, Object> DEFAULT_CLASSIFIER_PARAMETERS = new HashMap<String, Object>() {
        {
            put(ClassifierTuner.NUM_TRAINING_ROUNDS_KEY, DEFAULT_NUM_TRAINING_ROUNDS);

            put(LEARNING_RATE_KEY, DEFAULT_LEARNING_RATE);
            put(MAX_DEPTH_KEY, DEFAULT_MAX_DEPTH);

            put(OBJECTIVE_KEY, DEFAULT_OBJECTIVE_VALUE);
            put(EVAL_METRIC_KEY, DEFAULT_EVAL_METRIC);

            put(SILENT_KEY, DEFAULT_SILENT);
            put(SEED_KEY, DEFAULT_SEED);
        }
    };
    @SuppressWarnings("serial")
    public static final Map<String, ClassifierParamRange<?>> DEFAULT_TUNING_PARAMETERS
            = new HashMap<String, ClassifierParamRange<?>>() {
        {
            put(LEARNING_RATE_KEY, new ClassifierParamRange.ClassifierLogParamRange(0.01, 10.0));
            put(MAX_DEPTH_KEY, new ClassifierParamRange.ClassifierIntegerLinearParamRange(2, 20));
            put(GAMMA_KEY, new ClassifierParamRange.ClassifierLinearParamRange(0.0, 20.0));
            put(MIN_CHILD_WEIGHT_KEY, new ClassifierParamRange.ClassifierLogParamRange(1.0, 100.0));
            put(SUBSAMPLE_KEY, new ClassifierParamRange.ClassifierLinearParamRange(0.5, 1.0));
            put(COLSAMPLE_BY_TREE_KEY, new ClassifierParamRange.ClassifierLinearParamRange(0.5, 1.0));
            put(COLSAMPLE_BY_LEVEL_KEY, new ClassifierParamRange.ClassifierLinearParamRange(0.5, 1.0));
            put(MAX_DELTA_STEP_KEY, new ClassifierParamRange.ClassifierLinearParamRange(0.0, 10.0));
        }
    };
    @SuppressWarnings("serial")
    public static final Map<String, Boolean> BUILTIN_MAXIMIZE_EVAL_METRIC = new HashMap<String, Boolean>() {
        {
            put("rmse", false);
            put("mae", false);
            put("logloss", false);
            put("error", false);
            put("merror", false);
            put("auc", true);
            put("ndcg", true);
            put("ndcg-", true);
            put("map", true);
            put("map-", true);
            put("poisson-nloglik", false);
            put("gamma-nloglik", false);
            put("cox-nloglik", false);
            put("gamma-deviance", false);
            put("tweedie-nloglik", false);
        }
    };

}
