package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.AlignmentWrapperFactory
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class RecordAlignmentWrapperFactory implements AlignmentWrapperFactory<Record, RecordAlignment> {
    @Override
    RecordAlignment create(Record tcrData1, Record tcrData2, Alignment<AminoAcidSequence> alignment) {
        new RecordAlignment(tcrData1, tcrData2, alignment)
    }
}
