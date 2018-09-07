package org.broadinstitute.hellbender.tools.spark.sv.ml;

import htsjdk.samtools.util.IOUtil;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple class to hold truth sets: matrix of features paired with array of class labels
 * Note: if classLabels is passed in as null, it will be converted to an empty array so that checks for presence of
 * classLabels can just check .length
 * */
public class TruthSet {
    public static final int CLASS_LABEL_COLUMN = 0;

    private final RealMatrix features;
    private final int[] classLabels;

    TruthSet(final RealMatrix features, final int[] classLabels) {
        if(classLabels != null && classLabels.length > 0 && features.getRowDimension() != classLabels.length) {
            throw new IllegalArgumentException(
                    "classLabels must either be null, empty, or have same length as the number of rows in features."
            );
        }
        this.features = features;
        this.classLabels = classLabels == null ? new int[0] : classLabels;
    }

    /**
     * Load TruthSet from (possibly gzipped) csv file. The file must represent numeric types only.
     * Lines beginning with "#" will be ignored.
     * Values in CLASS_LABEL_COLUMN of file will be interpreted as class labels.
     */
    public static TruthSet loadCsvFile(final String filename) {
        return loadCsvFile(filename, ",", "#", CLASS_LABEL_COLUMN);
    }

    /**
     * Load TruthSet from (possibly gzipped) file which uses arbitrary delimiter to sparate values.
     * @param delimiter specify separator for values in file
     * @param commentCharacter lines beginning with this value are ignored
     * @param classLabelsColumn values in this column of the file will be interpreted as class labels.
     *                          If classLabelsColumn < 0, then classLabels will be empty
     */
    public static TruthSet loadCsvFile(final String filename, final String delimiter,
                                       final String commentCharacter, final int classLabelsColumn) {
        final int numColumns;
        final List<double[]> rowsList = new ArrayList<>();
        final List<Integer> classLabelsList = new ArrayList<>();
        try (final BufferedReader reader = IOUtil.openFileForBufferedReading(new File(filename))) {
            if(!reader.ready()) {
                throw new GATKException("Unable to read matrix from " + filename);
            }
            getNextCsvFeaturesRow(reader, delimiter, commentCharacter, classLabelsColumn, rowsList, classLabelsList, -1);
            numColumns = rowsList.get(0).length;
            while(reader.ready()) {
                getNextCsvFeaturesRow(reader, delimiter, commentCharacter, classLabelsColumn, rowsList, classLabelsList, numColumns);
            }
        } catch(IOException err) {
            throw new GATKException(err.getMessage());
        }
        return new TruthSet(
                new Array2DRowRealMatrix(rowsList.toArray(new double[0][]), false),
                classLabelsList.stream().mapToInt(Integer::intValue).toArray()
        );
    }

    private static void getNextCsvFeaturesRow(final BufferedReader reader, final String delimiter,
                                              final String commentCharacter, final int classLabelsColumn,
                                              final List<double[]> rowsList, final List<Integer> classLabelsList,
                                              final int numColumns) throws IOException {
        String line = reader.readLine();
        while(line.startsWith(commentCharacter) || line.isEmpty()) {
            line = reader.readLine();
        }
        final String[] words = line.split(delimiter, -1);
        final int numColumnsInRow = classLabelsColumn < 0 ? words.length : words.length - 1;
        if(numColumns >= 0 && numColumnsInRow != numColumns) {
            throw new GATKException("filename does not encode a matrix, rows do not all have the same length");
        }
        final double[] features = new double[numColumnsInRow];
        int featureIndex = 0;
        for(int columnNumber = 0; columnNumber < words.length; ++columnNumber) {
            final String word = words[columnNumber];
            if(columnNumber == classLabelsColumn) {
                classLabelsList.add(Integer.valueOf(word));
            } else {
                features[featureIndex] = Double.valueOf(word);
                ++featureIndex;
            }
        }
        rowsList.add(features);
    }

    public int getNumRows() {
        return features.getRowDimension();
    }

    public int getNumColumns() {
        return features.getColumnDimension();
    }

    public TruthSet sliceRows(final int[] rowIndices) {
        return new TruthSet(
                MachineLearningUtils.sliceRows(features, rowIndices),
                MachineLearningUtils.slice(classLabels, rowIndices)
        );
    }

    /**
     * Extract a stratify array from a TruthSet. Columns are binned and stratify array values are obtained by finding
     * unique rows.
     * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
     * @param numBins Number of bins into which to group values in each column
     * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
     *                                  the number of columns in consideration is reduced and the under-count rows are
     *                                  searched again for unique values.
     * @return stratifyArray, int array with length = number of rows in features
     */
    public int[] getStratifyArray(final int numBins, final int minCountsPerStratifyValue) {
        return getStratifyArray(numBins, minCountsPerStratifyValue, new HashSet<>());
    }

    /**
     * Extract a stratify array from a TruthSet. Columns are binned and stratify array values are obtained by finding
     * unique rows. This version can specify some columns as categorical (unbinnable).
     * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
     * @param numBins Number of bins into which to group values in each column
     * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
     *                                  the number of columns in consideration is reduced and the under-count rows are
     *                                  searched again for unique values.
     * @param categoricalColumns Columns in categoricalColumns are not binned.
     * @return stratifyArray, int array with length = number of rows in features
     */
    public int[] getStratifyArray(final int numBins, final int minCountsPerStratifyValue,
                                  final Set<Integer> categoricalColumns) {
        // class labels are the most important stratify column, so make them first. Also make them categorical
        final RealMatrix stratifyMatrix = MachineLearningUtils.concatenateColumns(
                new Array2DRowRealMatrix(IntStream.of(classLabels).mapToDouble(Double::valueOf).toArray()),
                features
        );
        final Set<Integer> concatenatedCategoricalColumns = categoricalColumns == null ? new HashSet<>()
                : categoricalColumns.stream().map(x -> x + 1).collect(Collectors.toSet());
        concatenatedCategoricalColumns.add(0); // column 0 (classLabel) is categorical
        return MachineLearningUtils.stratifyMatrixToStratifyArray(stratifyMatrix, numBins, minCountsPerStratifyValue, concatenatedCategoricalColumns);
    }

    public RealMatrix getFeatures() {
        return features;
    }

    public int[] getClassLabels() {
        return classLabels;
    }
}
