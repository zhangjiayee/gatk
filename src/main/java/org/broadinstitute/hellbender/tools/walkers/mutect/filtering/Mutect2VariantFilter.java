package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.IndexRange;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class Mutect2VariantFilter extends Mutect2Filter {
    public Mutect2VariantFilter() { }

    @Override
    public List<Double> errorProbabilities(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        int numAltAlleles = vc.getNAlleles() - 1;
        final double result = Mutect2FilteringEngine.roundFinitePrecisionErrors(requiredAnnotations().stream().allMatch(vc::hasAttribute) ?
                calculateErrorProbability(vc, filteringEngine, referenceContext) : 0.0);
        return Collections.nCopies(numAltAlleles, result);
    }

    protected abstract double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);
}
