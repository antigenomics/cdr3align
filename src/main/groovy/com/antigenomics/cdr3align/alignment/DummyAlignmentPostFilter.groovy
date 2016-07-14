package com.antigenomics.cdr3align.alignment

import com.milaboratory.core.alignment.Alignment

class DummyAlignmentPostFilter implements AlignmentPostFilter {
    @Override
    boolean pass(Alignment alignment) {
        true
    }
}
