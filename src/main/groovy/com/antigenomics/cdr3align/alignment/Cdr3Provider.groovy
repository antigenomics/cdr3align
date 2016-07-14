package com.antigenomics.cdr3align.alignment

import com.milaboratory.core.sequence.AminoAcidSequence

interface Cdr3Provider<T> {
    AminoAcidSequence getCdr3Aa(T tcrData)
}