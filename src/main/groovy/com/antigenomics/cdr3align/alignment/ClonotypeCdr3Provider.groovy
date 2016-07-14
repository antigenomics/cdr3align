package com.antigenomics.cdr3align.alignment

import com.antigenomics.vdjtools.sample.Clonotype
import com.milaboratory.core.sequence.AminoAcidSequence

class ClonotypeCdr3Provider implements Cdr3Provider<Clonotype>{
    @Override
    AminoAcidSequence getCdr3Aa(Clonotype clonotype) {
        clonotype.cdr3aaBinary
    }
}
