package com.antigenomics.cdr3align.alignment

import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class DummyAlignmentWrapperFactory<T> implements AlignmentWrapperFactory<T, AlignmentWrapper<T>> {
    @Override
    AlignmentWrapper<T> create(T tcrData1, T tcrData2, Alignment<AminoAcidSequence> alignment) {
        new AlignmentWrapper<T>(tcrData1, tcrData2, alignment)
    }
}
