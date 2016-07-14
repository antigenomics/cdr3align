package com.antigenomics.cdr3align.db

import com.antigenomics.cdr3align.alignment.AlignmentPostFilter
import com.milaboratory.core.alignment.Alignment

class DummyAlignmentPostFilter implements AlignmentPostFilter {
    @Override
    boolean pass(Alignment alignment) {
        true
    }
}
