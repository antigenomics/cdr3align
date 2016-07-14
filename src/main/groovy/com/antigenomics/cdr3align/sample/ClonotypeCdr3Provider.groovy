package com.antigenomics.cdr3align.sample

import com.antigenomics.cdr3align.alignment.Cdr3Provider
import com.antigenomics.vdjtools.sample.Clonotype
import com.milaboratory.core.sequence.AminoAcidSequence

class ClonotypeCdr3Provider implements Cdr3Provider<Clonotype> {
    @Override
    AminoAcidSequence getCdr3Aa(Clonotype clonotype) {
        clonotype.cdr3aaBinary
    }
}
