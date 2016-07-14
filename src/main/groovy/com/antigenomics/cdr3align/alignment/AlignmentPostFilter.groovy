package com.antigenomics.cdr3align.alignment

import com.milaboratory.core.alignment.Alignment

interface AlignmentPostFilter {
    boolean pass(Alignment alignment)
}