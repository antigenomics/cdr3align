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

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Record record = (Record) o

        if (cdr3 != record.cdr3) return false

        return true
    }

    int hashCode() {
        return cdr3.hashCode()
    }
}
