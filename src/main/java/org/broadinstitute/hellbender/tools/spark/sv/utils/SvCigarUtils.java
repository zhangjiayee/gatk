package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
public final class SvCigarUtils {

    /**
     * Computes the corresponding distance needs to be walked on the read, given the Cigar and distance walked on the reference.
     * @param cigar   cigar along the 5-3 direction of read (when read is mapped to reverse strand, bwa mem output cigar should be inverted)
     * @param start      start position (1-based) on the read (note it should not count the hard clipped bases, as usual)
     * @param refDist                 distance to walk on the reference
     * @param backward              whether to walk backwards along the read or not
     * @return                          corresponding walk distance on read (always positive)
     * @throws IllegalArgumentException if input cigar contains padding operation or 'N', or
     *                                  either of {@code start} or distance is non-positive, or
     *                                  {@code start} is larger than read length, or
     *                                  requested reference walk distance is longer than the total read bases in cigar, or
     *                                  computed read walk distance would "walk off" the read
     */
    @VisibleForTesting
    public static int computeAssociatedDistOnRead(final Cigar cigar, final int start, final int refDist, final boolean backward) {

        Utils.validateArg(refDist > 0 && start > 0, () -> "start " + start + " or distance " + refDist + " is non-positive.");

        final List<CigarElement> elements = backward ? Lists.reverse(cigar.getCigarElements()) : cigar.getCigarElements();

        final int readLength = elements.stream().mapToInt(ce -> ce.getOperator().consumesReadBases() ? ce.getLength() : 0).sum();
        final int readBasesToSkip = backward ? readLength - start : start - 1;

        int readBasesConsumed = 0;
        int refBasesConsumed = 0;

        for (final CigarElement element : elements){
            final int readBasesConsumedBeforeElement = readBasesConsumed;

            readBasesConsumed += element.getOperator().consumesReadBases() ? element.getLength() : 0;
            // skip cigar elements that end before the read start or start after the reference end
            if (readBasesConsumed <= readBasesToSkip) {
                continue;
            }

            refBasesConsumed += element.getOperator().consumesReferenceBases() ? element.getLength() - Math.max(readBasesToSkip - readBasesConsumedBeforeElement, 0) : 0;
            if (refBasesConsumed >= refDist) {
                final int excessRefBasesInElement = Math.max(refBasesConsumed - refDist, 0);
                return readBasesConsumed - readBasesToSkip - (element.getOperator().consumesReadBases() ? excessRefBasesInElement : 0);
            }
        }

        throw new IllegalArgumentException("Cigar " + cigar + "does not contain at least " + refDist + " reference bases past red start " + start + ".");
    }
}
