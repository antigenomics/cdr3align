package com.antigenomics.cdr3align.db

import com.milaboratory.core.sequence.AminoAcidSequence

class Record {
    final AminoAcidSequence cdr3
    final List<Segment> vSegments, jSegments
    final Gene gene
    final List<AminoAcidSequence> antigen = new ArrayList<>()

    Record(String gene, String cdr3, String vSegments, String jSegments) {
        this.gene = Gene."$gene"
        this.cdr3 = new AminoAcidSequence(cdr3)
        this.vSegments = vSegments.split(",").collect { SegmentCache.INSTANCE.getOrCreate(it) }
        this.jSegments = jSegments.split(",").collect { SegmentCache.INSTANCE.getOrCreate(it) }
    }

    @Override
    boolean equals(o) {
        this.cdr3 == ((Record) o).cdr3
    }

    @Override
    int hashCode() {
        cdr3.hashCode()
    }
}
