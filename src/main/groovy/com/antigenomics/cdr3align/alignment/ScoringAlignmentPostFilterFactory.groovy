package com.antigenomics.cdr3align.alignment

import com.antigenomics.vdjdb.scoring.AlignmentScoring


class ScoringAlignmentPostFilterFactory<T> implements AlignmentPostFilterFactory<T> {
    final AlignmentScoring scoring

    ScoringAlignmentPostFilterFactory(AlignmentScoring scoring) {
        this.scoring = scoring
    }

    @Override
    AlignmentPostFilter getPostFilter(Cdr3Provider<T> cdr3Provider, T tcrData) {
        def cdr3 = cdr3Provider.getCdr3Aa(tcrData)
        new ScoringAlignmentPostFilter(scoring.computeBaseScore(cdr3), cdr3.size(), scoring)
    }
}
