package com.antigenomics.cdr3align.alignment

interface AlignmentPostFilterFactory<T> {
    AlignmentPostFilter getPostFilter(Cdr3Provider<T> cdr3Provider, T tcrData)
}