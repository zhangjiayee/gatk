package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
public final class SvCigarUtils {

    /**
     * Checks the input CIGAR for assumption that operator 'D' is not immediately adjacent to clipping operators.
     * Then convert the 'I' CigarElement, if it is at either end (terminal) of the input cigar, to a corresponding 'S' operator.
     * Note that we allow CIGAR of the format '10H10S10I10M', but disallows the format if after the conversion the cigar turns into a giant clip,
     * e.g. '10H10S10I10S10H' is not allowed (if allowed, it becomes a giant clip of '10H30S10H' which is non-sense).
     *
     * @return a pair of number of clipped (hard and soft, including the ones from the converted terminal 'I') bases at the front and back of the
     *         input {@code cigarAlongInput5to3Direction}.
     *
     * @throws IllegalArgumentException when the checks as described above fail.
     */
    @VisibleForTesting
    public static Cigar checkCigarAndConvertTerminalInsertionToSoftClip(final Cigar cigar) {

        if (cigar.numCigarElements()<2 ) return cigar.getCigarElements();

        final List<CigarElement> cigarElements = new ArrayList<>(cigar.getCigarElements());

        final List<CigarElement> convertedList = convertInsToSoftClipFromOneEnd(cigarElements, true);
        return convertInsToSoftClipFromOneEnd(convertedList, false);
    }

    /**
     * Computes the corresponding distance needs to be walked on the reference, given the Cigar and distance walked on the read.
     * @param cigarAlong5To3DirOfRead   cigar along the 5-3 direction of read (when read is mapped to reverse strand, bwa mem output cigar should be inverted)
     * @param startInclusiveOnRead      start position (1-based) on the read (note it should not count the hard clipped bases, as usual)
     * @param distanceOnRead            distance to walk on the read
     * @return                          corresponding walk distance on reference
     * @throws IllegalArgumentException if input cigar contains padding operation or 'N', or
     *                                  either of startInclusive or distance is non-positive, or
     *                                  startInclusive + distance -1 is longer than the read as suggested by the cigar
     */
    @VisibleForTesting
    public static int computeAssociatedDistOnRef(final Cigar cigarAlong5To3DirOfRead, final int startInclusiveOnRead,
                                                 final int distanceOnRead) {

        final int endInclusive = startInclusiveOnRead + distanceOnRead - 1;
        Utils.validateArg(startInclusiveOnRead>0 && distanceOnRead > 0,
                "start position (" + startInclusiveOnRead + ") or distance (" + distanceOnRead + ") is non-positive.");
        final List<CigarElement> cigarElements = cigarAlong5To3DirOfRead.getCigarElements();
        Utils.validateArg(cigarElements.stream().noneMatch(ce -> ce.getOperator().isPadding() || ce.getOperator().equals(CigarOperator.N)),
                "cigar contains padding, which is currently unsupported; cigar: " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));
        Utils.validateArg(cigarElements.stream()
                        .mapToInt(ce -> ce.getOperator().consumesReadBases() ? ce.getLength() : 0).sum() >= endInclusive,
                "start location (" + startInclusiveOnRead + ") and walking distance (" + distanceOnRead +
                        ") would walk out of the read, indicated by cigar " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));

        int readBasesConsumed = 0;

        // skip first several elements that give accumulated readBasesConsumed below startInclusiveOnRead
        int idx = 0;
        CigarElement currEle = cigarElements.get(idx);
        while (readBasesConsumed + (currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0) < startInclusiveOnRead) {
            readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
            currEle = cigarElements.get(++idx);
        }
        int refWalkDist = 0;
        int readWalked = 0;
        while (idx != cigarElements.size()) {
            currEle = cigarElements.get(idx);
            final int skip = Math.max(0, startInclusiveOnRead - readBasesConsumed - 1);
            final int effectiveLen = currEle.getLength() - skip;

            if (currEle.getOperator().consumesReadBases()) {
                if (readWalked + effectiveLen < distanceOnRead) { // hasn't walked enough yet on read
                    readWalked += effectiveLen;
                    refWalkDist += currEle.getOperator().consumesReferenceBases() ? effectiveLen : 0;
                    readBasesConsumed += currEle.getOperator().consumesReferenceBases() ? currEle.getLength() : 0;
                } else { // would be walked enough on read
                    refWalkDist += currEle.getOperator().consumesReferenceBases() ? distanceOnRead - readWalked : 0;
                    readWalked = distanceOnRead;
                    break;
                }
            } else {
                refWalkDist += currEle.getOperator().consumesReferenceBases() ? effectiveLen : 0;
                readBasesConsumed += currEle.getOperator().consumesReferenceBases() ? currEle.getLength() : 0;
            }

            ++idx;
        }

        return refWalkDist;
    }

    /**
     * Computes the corresponding distance needs to be walked on the read, given the Cigar and distance walked on the reference.
     * @param cigarAlong5To3DirOfRead   cigar along the 5-3 direction of read (when read is mapped to reverse strand, bwa mem output cigar should be inverted)
     * @param startInclusiveOnRead      start position (1-based) on the read (note it should not count the hard clipped bases, as usual)
     * @param distOnRef                 distance to walk on the reference
     * @param walkBackward              whether to walk backwards along the read or not
     * @return                          corresponding walk distance on read (always positive)
     * @throws IllegalArgumentException if input cigar contains padding operation or 'N', or
     *                                  either of {@code startInclusiveOnRead} or distance is non-positive, or
     *                                  {@code startInclusiveOnRead} is larger than read length, or
     *                                  requested reference walk distance is longer than the total read bases in cigar, or
     *                                  computed read walk distance would "walk off" the read
     */
    @VisibleForTesting
    public static int computeAssociatedDistOnRead(final Cigar cigarAlong5To3DirOfRead, final int startInclusiveOnRead,
                                                  final int distOnRef, final boolean walkBackward) {

        Utils.validateArg(distOnRef > 0 && startInclusiveOnRead > 0,
                "start position (" + startInclusiveOnRead + ") or distance (" + distOnRef + ") is non-positive.");

        final List<CigarElement> cigarElementsInOrderOfWalkingDir = (walkBackward ? CigarUtils.invertCigar(cigarAlong5To3DirOfRead): cigarAlong5To3DirOfRead).getCigarElements();
        Utils.validateArg(cigarElementsInOrderOfWalkingDir.stream().noneMatch(ce -> ce.getOperator().isPadding() || ce.getOperator().equals(CigarOperator.N)),
                "cigar contains padding, which is currently unsupported; cigar: " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));
        final int readUnclippedLength = CigarUtils.countUnclippedReadBases(cigarAlong5To3DirOfRead);
        Utils.validateArg(readUnclippedLength >= startInclusiveOnRead,
                "given start location on read (" + startInclusiveOnRead + ") is higher than read unclipped length (" + readUnclippedLength+ "), cigar: " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));
        final int totalRefLen = cigarElementsInOrderOfWalkingDir.stream().mapToInt(ce -> ce.getOperator().consumesReferenceBases() ? ce.getLength() : 0).sum();
        Utils.validateArg(totalRefLen >= distOnRef,
                "given walking distance on reference (" + distOnRef + ") would be longer than the total number (" +
                        + totalRefLen + ") of reference bases spanned by the cigar, indicated by cigar " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));

        final int readLength = cigarElementsInOrderOfWalkingDir.stream().mapToInt(ce -> ce.getOperator().consumesReadBases() ? ce.getLength() : 0).sum();
        final int effectiveReadStartInclusive = walkBackward ? readLength - startInclusiveOnRead + 1 : startInclusiveOnRead;
        // skip first several elements that give accumulated readBasesConsumed below startInclusiveOnRead
        int idx = 0;
        int readBasesConsumed = 0;
        CigarElement currEle = cigarElementsInOrderOfWalkingDir.get(idx);
        while (readBasesConsumed + (currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0) < effectiveReadStartInclusive) {
            readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
            currEle = cigarElementsInOrderOfWalkingDir.get(++idx);
        }

        // when we reach here, we have skipped just enough read bases to start counting ref bases or, currEle would lead us to such state
        int readWalkDist = 0;
        int refWalked = 0;
        while (idx != cigarElementsInOrderOfWalkingDir.size()) {
            currEle = cigarElementsInOrderOfWalkingDir.get(idx);
            final int skip = Math.max(0, effectiveReadStartInclusive - readBasesConsumed - 1);
            final int effectiveLen = currEle.getLength() - skip;

            if (currEle.getOperator().consumesReferenceBases()) {
                if (refWalked + effectiveLen < distOnRef) { // hasn't walked enough yet on reference
                    refWalked += effectiveLen;
                    readWalkDist += currEle.getOperator().consumesReadBases() ? effectiveLen : 0;
                    readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
                } else { // would be walked enough on reference
                    readWalkDist += currEle.getOperator().consumesReadBases() ? distOnRef - refWalked : 0;
                    refWalked = distOnRef;
                    break;
                }
            } else {
                readWalkDist += currEle.getOperator().consumesReadBases() ? effectiveLen : 0;
                readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
            }
            ++idx;
        }

        if (refWalked < distOnRef)
            throw new IllegalArgumentException("Computed walk distance (start: " + startInclusiveOnRead + ", distOnRef: " +
                    distOnRef + ") on read beyond read length (" + readLength +") with cigar " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));

        return readWalkDist;
    }

}
