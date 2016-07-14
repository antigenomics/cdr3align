package com.antigenomics.cdr3align.alignment

class DummyAlignmentPostFilterFactory<T> implements AlignmentPostFilterFactory<T> {
    static final DummyAlignmentPostFilter DUMMY_FILTER_INSTANCE = new DummyAlignmentPostFilter()

    @Override
    AlignmentPostFilter getPostFilter(Cdr3Provider<T> cdr3Provider, T tcrData) {
        DUMMY_FILTER_INSTANCE
    }
}
