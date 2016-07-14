package com.antigenomics.cdr3align.sample

import com.antigenomics.cdr3align.alignment.AlignmentPostFilter
import com.antigenomics.cdr3align.alignment.AlignmentPostFilterFactory
import com.antigenomics.cdr3align.alignment.Cdr3Provider
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
