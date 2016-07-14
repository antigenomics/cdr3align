package com.antigenomics.cdr3align.sample

import com.antigenomics.cdr3align.alignment.AlignmentWrapper
import com.antigenomics.cdr3align.alignment.BatchAligner
import com.antigenomics.vdjdb.scoring.AlignmentScoring
import com.antigenomics.vdjtools.sample.Clonotype

class SampleBatchAligner extends BatchAligner<Clonotype, AlignmentWrapper<Clonotype>> {
    SampleBatchAligner(int mismatches, int deletions, int insertions, int total,
                       AlignmentScoring alignmentScoring) {
        super(mismatches, deletions, insertions, total,
                new ClonotypeCdr3Provider(),
                new DummyAlignmentWrapperFactory<>(),
                new DummyAlignmentPreFilter<>(),
                new ScoringAlignmentPostFilterFactory<>(alignmentScoring))
    }
}
