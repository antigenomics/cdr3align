package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.BatchAligner

class VdjdbBatchAligner extends BatchAligner<Record, RecordAlignment> {
    VdjdbBatchAligner(int mismatches, int deletions, int insertions, int total,
                      boolean requireVMatch, boolean requireJMatch) {
        super(mismatches, deletions, insertions, total,
                new RecordCdr3Provider<>(),
                new RecordAlignmentWrapperFactory<>(),
                new RecordAlignmentPreFilter<>(requireVMatch, requireJMatch),
                new DummyAlignmentPostFilterFactory<>())
    }
}
