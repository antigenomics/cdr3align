package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.AlignmentPreFilter

class RecordAlignmentPreFilter implements AlignmentPreFilter<Record> {
    @Override
    boolean canAlign(Record tcrData1, Record tcrData2) {
        tcrData1.gene == tcrData2.gene
    }
}
