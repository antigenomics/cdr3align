package com.antigenomics.cdr3align.alignment

interface AlignmentPreFilter<T> {
    boolean canAlign(T tcrData1, T tcrData2)
}
