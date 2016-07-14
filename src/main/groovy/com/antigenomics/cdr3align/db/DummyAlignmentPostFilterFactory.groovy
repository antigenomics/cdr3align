package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.AlignmentPostFilter
import com.antigenomics.cdr3align.alignment.AlignmentPostFilterFactory
import com.antigenomics.cdr3align.alignment.Cdr3Provider

class DummyAlignmentPostFilterFactory<T> implements AlignmentPostFilterFactory<T> {
    static final DummyAlignmentPostFilter DUMMY_FILTER_INSTANCE = new DummyAlignmentPostFilter()

    @Override
    AlignmentPostFilter getPostFilter(Cdr3Provider<T> cdr3Provider, T tcrData) {
        DUMMY_FILTER_INSTANCE
    }
}
