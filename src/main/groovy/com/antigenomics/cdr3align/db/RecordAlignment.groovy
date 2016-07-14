package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.AlignmentWrapper
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class RecordAlignment extends AlignmentWrapper<Record> {
    final boolean antigensMatch

    RecordAlignment(Record record1, Record record2, Alignment<AminoAcidSequence> alignment) {
        super(record1, record2, alignment)
        this.antigensMatch = record1.antigen.any { record2.antigen.contains(it) }
    }
}
