package com.antigenomics.cdr3align.alignment

import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

interface AlignmentWrapperFactory<T, R extends AlignmentWrapper<T>> {
    R create(T tcrData1, T tcrData2, Alignment<AminoAcidSequence> alignment)
}