package org.broadinstitute.hellbender.tools.spark.sv.ml;

import org.apache.commons.lang3.StringUtils;

import java.util.Random;

/**
 * Interface for classifier parameter range.
 */
public interface ClassifierParamRange<T> {
    T[] getRandomSamples(final Random random, final int numSamples);

    /**
     * Class for double-valued ClassifierParamRange with linear (uniform) sampling across range
     */
    class ClassifierLinearParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        public ClassifierLinearParamRange(double low, double high) {
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Double[] getRandomSamples(final Random random, final int numSamples) {
            final Double[] samples = new Double[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (high + low) / 2.0;
                }
                return samples;
            }
            final double delta = (high - low) / (numSamples - 1);
            double val = low;
            final int[] permutation = MachineLearningUtils.getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val += delta;
                samples[permutation[i]] = val;
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Class for double-valued ClassifierParamRange with logarithmic sampling across range
     */
    class ClassifierLogParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        public ClassifierLogParamRange(double low, double high) {
            if(low * high <= 0) {
                throw new IllegalArgumentException("low (" + low + ") and high (" + high + ") must be the same sign, and non-zero");
            }
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Double[] getRandomSamples(final Random random, final int numSamples) {
            final Double[] samples = new Double[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = Math.sqrt(high * low);
                }
                return samples;
            }
            final double delta = Math.pow(high / low, 1.0 / (numSamples - 1));
            double val = low;
            final int[] permutation = MachineLearningUtils.getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val *= delta;
                samples[permutation[i]] = val;
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Class for integer-valued ClassifierParamRange with linear (uniform) sampling across range
     */
    class ClassifierIntegerLinearParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        public ClassifierIntegerLinearParamRange(int low, int high) {
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Integer[] getRandomSamples(final Random random, final int numSamples) {
            final Integer[] samples = new Integer[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (int)Math.round((high + low) / 2.0);
                }
                return samples;
            }
            final double delta = (high - low) / (double)(numSamples - 1);
            double val = low;
            final int[] permutation = MachineLearningUtils.getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val += delta;
                samples[permutation[i]] = (int)Math.round(val);
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Class for integer-valued ClassifierParamRange with logarithmic sampling across range
     */
    class ClassifierIntegerLogParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        public ClassifierIntegerLogParamRange(int low, int high) {
            if(low * high <= 0) {
                throw new IllegalArgumentException("low (" + low + ") and high (" + high + ") must be the same sign, and non-zero");
            }
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Integer[] getRandomSamples(final Random random, final int numSamples) {
            final Integer[] samples = new Integer[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (int)Math.round(Math.sqrt((double)high * low));
                }
                return samples;
            }
            final double delta = Math.pow(high / (double)low, 1.0 / (numSamples - 1));
            double val = low;
            final int[] permutation = MachineLearningUtils.getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val *= delta;
                samples[permutation[i]] = (int)Math.round(val);
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Progress bar for console programs; to show work done, work to-do, elapsed time, and estimated remaining time for
     * long-running tasks. It avoids redrawing too frequently, so as to prevent slowing progress of actual work.
     */
    class ConsoleProgressBar {
        private static final long MIN_UPDATE_INTERVAL_NS = 500000000; // = 0.5 sec
        private static final int NUM_BAR_CHARACTERS = 10;

        private static final double SECONDS_IN_MINUTE = 60.0;
        private static final long MINUTES_IN_HOUR = 60;
        private static final long HOURS_IN_DAY = 24;
        private static final String CARRIAGE_RETURN = "\r";

        private final long workToDo;
        private final long bornTime;
        private long nextUpdateTime;
        private long workDone;
        private long workRemaining;
        private int maxOutputLength;
        private final String updateInfoFormat;

        /**
         * Create progress bar
         * @param workToDo sum of amount of work for all tasks. Usually this will just be number of tasks.
         */
        ConsoleProgressBar(final long workToDo) {
            if(workToDo <= 0) {
                throw new IllegalArgumentException("workToDo must be > 0");
            }
            this.workToDo = workToDo;
            workDone = 0;
            workRemaining = workToDo;
            bornTime = System.nanoTime();
            nextUpdateTime = bornTime + MIN_UPDATE_INTERVAL_NS;
            maxOutputLength = 0;
            final int workToDoLength = String.format("%d", workToDo).length();
            final String workDoneFormat = String.format("%%%dd/%%%dd", workToDoLength, workToDoLength);
            updateInfoFormat = workDoneFormat + " %.1f%% elapsed %s, remaining %s";
            drawBar(0.0, Double.NaN);
        }

        /**
         * Update the progress bar. Re-draw if enough time has elapsed.
         * @param workJustCompleted This value can be altered if different tasks have quantifiable differences in the
         *                          amount of work (or time) they will take. Probably this will normally just equal 1.
         */
        public void update(final long workJustCompleted) {
            if(workJustCompleted <= 0) {
                throw new IllegalArgumentException("workJustCompleted must be > 0");
            }
            workRemaining -= workJustCompleted;
            workDone += workJustCompleted;
            final long now = System.nanoTime();
            if(now < nextUpdateTime && workRemaining > 0) {
                return;  // avoid thrashing to the screen
            } else {
                nextUpdateTime = now + MIN_UPDATE_INTERVAL_NS;
            }
            final double elapsedTimeSec = 1.0e-9 * (now - bornTime);
            final double workPerSec = workDone / elapsedTimeSec;
            final double remainingTimeSec = workRemaining / workPerSec;
            System.out.flush();
            drawBar(elapsedTimeSec, remainingTimeSec);
        }

        private void drawBar(final double elapsedTimeSec, final double remainingTimeSec) {
            // NOTE: carriage return means each output will start from the beginning of the line
            final StringBuilder stringBuilder = new StringBuilder(CARRIAGE_RETURN);
            // draw actual progress bar
            final int numBarFilled = (int)(NUM_BAR_CHARACTERS * workDone / workToDo);
            final int numBarUnfilled = NUM_BAR_CHARACTERS - numBarFilled;
            stringBuilder.append("|");
            stringBuilder.append(StringUtils.repeat('#', numBarFilled));
            stringBuilder.append(StringUtils.repeat(' ', numBarUnfilled));
            stringBuilder.append("| ");
            // write summary statistics on completion amount, times
            stringBuilder.append(
                    workDone > 0 ?
                    String.format(
                        updateInfoFormat, workDone, workToDo, workDone * 100.0 / workToDo,
                        secondsToTimeString(elapsedTimeSec), secondsToTimeString(remainingTimeSec)
                    )
                    : String.format(updateInfoFormat, 0, workToDo, 0.0, secondsToTimeString(0.0), "???")
            );
            // do any necessary padding
            if(stringBuilder.length() > maxOutputLength) {
                maxOutputLength = stringBuilder.length();
            } else if(stringBuilder.length() < maxOutputLength) {
                // pad with spaces to obliterate previous message
                stringBuilder.append(StringUtils.repeat(' ', maxOutputLength - stringBuilder.length()));
            }
            // if we're done, add newline
            if(workRemaining <= 0) {
                stringBuilder.append("\n");
            }
            // write out and flush
            System.out.print(stringBuilder.toString());
            System.out.flush();
        }

        private static String secondsToTimeString(double seconds) {
            if(seconds < SECONDS_IN_MINUTE) {
                return String.format("%.1fs", seconds);
            }
            long minutes = (int)Math.floor(seconds / SECONDS_IN_MINUTE);
            seconds = seconds % SECONDS_IN_MINUTE;
            long hours = minutes / MINUTES_IN_HOUR;
            if(hours <= 0) {
                return String.format("%dm %.1fs", minutes, seconds);
            }
            minutes = minutes % MINUTES_IN_HOUR;
            long days = hours / HOURS_IN_DAY;
            if(days <= 0) {
                return String.format("%dh %dm %.1fs", hours, minutes, seconds);
            } else {
                hours = hours % HOURS_IN_DAY;
                return String.format("%dd %dh %dm %.1fs", days, hours, minutes, seconds);
            }
        }
    }
}
