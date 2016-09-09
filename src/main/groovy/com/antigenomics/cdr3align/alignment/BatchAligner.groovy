package com.antigenomics.cdr3align.alignment

import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import groovyx.gpars.GParsPool

import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicInteger

class BatchAligner<T, R extends AlignmentWrapper<T>> {
    final TreeSearchParameters searchParameters
    final boolean allowSelfAlignments
    final Cdr3Provider<T> cdr3Provider
    final AlignmentWrapperFactory<T, R> alignmentWrapperFactory
    final AlignmentPreFilter<T> alignmentPreFilter
    final AlignmentPostFilterFactory<T> alignmentPostFilterFactory
    final SequenceTreeMap<AminoAcidSequence, T> treeMap = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET)

    BatchAligner(int mismatches, int deletions, int insertions, int total,
                 boolean allowSelfAlignments,
                 Cdr3Provider<T> cdr3Provider,
                 AlignmentWrapperFactory<T, R> alignmentWrapperFactory,
                 AlignmentPreFilter<T> alignmentPreFilter,
                 AlignmentPostFilterFactory<T> alignmentPostFilterFactory) {
        this.searchParameters = new TreeSearchParameters(mismatches, deletions, insertions, total)
        this.allowSelfAlignments = allowSelfAlignments
        this.cdr3Provider = cdr3Provider
        this.alignmentWrapperFactory = alignmentWrapperFactory
        this.alignmentPreFilter = alignmentPreFilter
        this.alignmentPostFilterFactory = alignmentPostFilterFactory
    }

    void add(Iterable<T> tcrData) {
        tcrData.each {
            treeMap.put(cdr3Provider.getCdr3Aa(it), it)
        }
    }

    Collection<R> align(int threads) {
        def alignments = new ConcurrentLinkedQueue<R>()
        def counter = new AtomicInteger()

        GParsPool.withPool(threads) {
            treeMap.values().findAll { alignmentPreFilter.notDummy(it) }.eachParallel { T from ->
                def fromCdr3 = cdr3Provider.getCdr3Aa(from)
                def iter = treeMap.getNeighborhoodIterator(fromCdr3, searchParameters)
                def to
                def alignmentPostFilter = alignmentPostFilterFactory.getPostFilter(cdr3Provider, from)

                def toCdr3Hash = new HashSet<AminoAcidSequence>()
                while ((to = iter.next()) != null) {
                    def toCdr3 = cdr3Provider.getCdr3Aa(to)
                    if ((allowSelfAlignments || fromCdr3 != toCdr3) &&
                            !toCdr3Hash.contains(toCdr3) &&
                            alignmentPreFilter.canAlign(from, to)) {
                        def alignment = iter.currentAlignment
                        if (alignmentPostFilter.pass(alignment)) {
                            alignments.add(alignmentWrapperFactory.create(from, to, alignment))
                            toCdr3Hash.add(toCdr3)
                        }
                    }
                }

                int counterVal
                if ((counterVal = counter.incrementAndGet()) % 100 == 0) {
                    println "[BatchAligner ${new Date().toString()}] Done all alignments for $counterVal records, " +
                            "total number of aligned CDR3 pairs is ${alignments.size()}"
                }
            }
        }

        alignments
    }
}
