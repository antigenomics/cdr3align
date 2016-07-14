package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.AlignmentPreFilter

class RecordAlignmentPreFilter implements AlignmentPreFilter<Record> {
    final boolean matchV, matchJ

    RecordAlignmentPreFilter(boolean matchV, boolean matchJ) {
        this.matchV = matchV
        this.matchJ = matchJ
    }

    @Override
    boolean canAlign(Record tcrData1, Record tcrData2) {
        tcrData1.gene == tcrData2.gene &&
                (!matchV || tcrData1.vSegments.any { tcrData2.vSegments.contains(it) }) &&
                (!matchJ || tcrData1.jSegments.any { tcrData2.jSegments.contains(it) })
    }
}
