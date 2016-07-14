package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.Cdr3Provider
import com.antigenomics.cdr3align.db.Record
import com.milaboratory.core.sequence.AminoAcidSequence

class RecordCdr3Provider implements Cdr3Provider<Record> {
    @Override
    AminoAcidSequence getCdr3Aa(Record record) {
        record.cdr3
    }
}
