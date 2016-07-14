package com.antigenomics.cdr3align.sample

import com.antigenomics.cdr3align.alignment.AlignmentPostFilter
import com.antigenomics.vdjdb.scoring.AlignmentScoring
import com.milaboratory.core.alignment.Alignment

class ScoringAlignmentPostFilter implements AlignmentPostFilter {
    final float baseScore
    final int refLength
    final AlignmentScoring scoring

    ScoringAlignmentPostFilter(float baseScore, int refLength,
                               AlignmentScoring scoring) {
        this.baseScore = baseScore
        this.refLength = refLength
        this.scoring = scoring
    }

    @Override
    boolean pass(Alignment alignment) {
        scoring.computeScore(alignment.absoluteMutations, baseScore, refLength)
    }
}
