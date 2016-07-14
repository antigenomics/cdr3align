package com.antigenomics.cdr3align.alignment

import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class AlignmentWrapper<T> {
    final T tcrData1, tcrData2
    final Alignment<AminoAcidSequence> alignment

    AlignmentWrapper(T tcrData1, T tcrData2, Alignment<AminoAcidSequence> alignment) {
        this.tcrData1 = tcrData1
        this.tcrData2 = tcrData2
        this.alignment = alignment
    }
}
