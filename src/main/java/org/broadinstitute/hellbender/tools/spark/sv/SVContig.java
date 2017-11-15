package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.ArraySVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.report.GATKReportColumnFormat;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

/**
 * Created by valentin on 9/16/17.
 */
@DefaultSerializer(SVContig.Serializer.class)
public class SVContig extends ArraySVHaplotype {

    private static final long serialVersionUID = 1L;

    public double getReferenceScore() {
        if (name.equals("ref")) {
            return 0.0;
        } else if (name.equals("alt")) {
            return Double.NEGATIVE_INFINITY;
        } else {
            return referenceAlignmentScore == null ? Double.NaN : referenceAlignmentScore.getValue();
        }
    }

    public double getAlternativeScore() {
        if (name.equals("alt")) {
            return 0.0;
        } else if (name.equals("ref")) {
            return Double.NEGATIVE_INFINITY;
        } else {
            return alternativeAlignmentScore == null ? Double.NaN  : alternativeAlignmentScore.getValue();
        }
    }

    enum Call {
        REF, ALT, NOCALL;
    }

    private final List<AlignmentInterval> referenceAlignment;
    private final List<AlignmentInterval> alternativeAlignment;
    private final AlignmentScore referenceAlignmentScore;
    private final AlignmentScore alternativeAlignmentScore;
    private final double callQuality;
    private final Call call;

    public static SVContig of(final GATKRead read) {
        final String variantId = getMandatoryAttribute(read, ComposeStructuralVariantHaplotypesSpark.VARIANT_CONTEXT_TAG);
        final List<AlignmentInterval> refAln = getAlignmentIntervalsAttribute(read, ComposeStructuralVariantHaplotypesSpark.REFERENCE_ALIGNMENT_TAG);
        final List<AlignmentInterval> altAln = getAlignmentIntervalsAttribute(read, ComposeStructuralVariantHaplotypesSpark.ALTERNATIVE_ALIGNMENT_TAG);
        final AlignmentScore refScore = getOptionalAlignmentScore(read, ComposeStructuralVariantHaplotypesSpark.REFERENCE_SCORE_TAG);
        final AlignmentScore altScore = getOptionalAlignmentScore(read, ComposeStructuralVariantHaplotypesSpark.ALTERNATIVE_SCORE_TAG);
        final SimpleInterval location = new SimpleInterval(read.getAssignedContig(), read.getAssignedStart(), read.getAssignedStart());
        return new SVContig(read.getName(), location, variantId, read.getBases(), refAln, refScore, altAln, altScore);

    }

    private static AlignmentScore getOptionalAlignmentScore(final GATKRead read, final String tag) {
        final String str = read.getAttributeAsString(tag);
        return str == null ? null : AlignmentScore.valueOf(str);
    }

    private static String getMandatoryAttribute(final GATKRead read, final String tag) {
        return ReadUtils.getOptionalStringAttribute(read, tag)
                .orElseThrow(() -> new UserException.BadInput("input read missing '" + tag + "' attribute"));
    }

    private static List<AlignmentInterval> getAlignmentIntervalsAttribute(final GATKRead read, final String tag) {
        final Optional<String> str = ReadUtils.getOptionalStringAttribute(read, tag);
        return str.isPresent() ? AlignmentInterval.decodeList(str.get()) : null;
    }

    public SVContig(final String name, final Locatable loc, final String variantId,
                     final byte[] bases, final List<AlignmentInterval> refAln, final AlignmentScore refScore,
                     final List<AlignmentInterval> altAln, final AlignmentScore altScore) {
        super(name, refAln, bases, variantId, SimpleInterval.of(loc), true);
        this.referenceAlignment = refAln;
        this.alternativeAlignment = altAln;
        this.referenceAlignmentScore = refScore;
        this.alternativeAlignmentScore = altScore;
        final double refScoreValue = refScore == null ? Double.NaN : refScore.getValue();
        final double altScoreValue = refScore == null ? Double.NaN : altScore.getValue();
        if (refScoreValue < altScoreValue) {
            this.call = Call.ALT;
            this.callQuality = altScoreValue - refScoreValue;
        } else if (altScoreValue < refScoreValue) {
            this.call = Call.REF;
            this.callQuality = refScoreValue - altScoreValue;
        } else {
            this.call = Call.NOCALL;
            this.callQuality = 0.0;
        }
    }

    private SVContig(final Kryo kryo, final Input input) {
        super(kryo, input);
        this.referenceAlignment = Serializer.readAlignment(kryo, input);
        this.alternativeAlignment = Serializer.readAlignment(kryo, input);
        this.referenceAlignmentScore = kryo.readObjectOrNull(input, AlignmentScore.class);
        this.alternativeAlignmentScore = kryo.readObjectOrNull(input, AlignmentScore.class);
        this.call = Call.valueOf(input.readString());
        this.callQuality = input.readDouble();
    }

    public static class Serializer extends ArraySVHaplotype.Serializer<SVContig> {

        @Override
        public void write(Kryo kryo, Output output, SVContig object) {
            super.write(kryo, output, object);
            writeAlignment(kryo, output, object.referenceAlignment);
            writeAlignment(kryo, output, object.alternativeAlignment);
            kryo.writeObjectOrNull(output, object.referenceAlignmentScore, AlignmentScore.class);
            kryo.writeObjectOrNull(output, object.alternativeAlignmentScore, AlignmentScore.class);
            output.writeString(object.call.name());
            output.writeDouble(object.callQuality);
        }

        @Override
        public SVContig read(Kryo kryo, Input input, Class<SVContig> type) {
            return new SVContig(kryo, input);
        }

        private void writeAlignment(final Kryo kryo, final Output output, final List<AlignmentInterval> alignment) {
            output.writeInt(alignment == null ? 0 : alignment.size());
            if (alignment != null) {
                for (final AlignmentInterval interval : alignment) {
                    kryo.writeObject(output, interval);
                }
            }
        }

        private static List<AlignmentInterval> readAlignment(final Kryo kryo, final Input input) {
            final int length = input.readInt();
            if (length <= 0) {
                return Collections.emptyList();
            } else {
                final AlignmentInterval[] intervals = new AlignmentInterval[length];
                for (int i = 0; i < length; i++) {
                    intervals[i] = kryo.readObject(input, AlignmentInterval.class);
                }
                return Collections.unmodifiableList(Arrays.asList(intervals));
            }
        }
    }

}