package com.antigenomics.cdr3align.sample

import com.antigenomics.cdr3align.alignment.AlignmentPreFilter

class DummyAlignmentPreFilter<T> implements AlignmentPreFilter<T> {
    @Override
    boolean canAlign(T tcrData1, T tcrData2) {
        true
    }

    @Override
    boolean notDummy(T tcrData) {
        true
    }
}
