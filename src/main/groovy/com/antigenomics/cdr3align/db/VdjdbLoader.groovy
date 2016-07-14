package com.antigenomics.cdr3align.db

import com.milaboratory.core.sequence.AminoAcidSequence

class VdjdbLoader {
    final String vdjdbSlimPath

    VdjdbLoader(String vdjdbSlimPath) {
        this.vdjdbSlimPath = vdjdbSlimPath
    }

    Collection<Record> load(int vdjdbConfThreshold,
                            List<String> genes, List<String> species,
                            int minCdr3PerAg) {
        def recordMap = new HashMap<String, Record>()

        def antigenCountMap = new HashMap<String, Integer>()

        def goodRecord = { List<String> record ->
            record[-1].toInteger() >= vdjdbConfThreshold &&
                    genes.any { it.equalsIgnoreCase(record[0]) } &&
                    species.any { it.equalsIgnoreCase(record[2]) }
        }

        def firstLine = true
        new File(vdjdbSlimPath).splitEachLine("\t") {
            if (!firstLine) {
                if (goodRecord(it)) {
                    antigenCountMap.put(it[3], (antigenCountMap[it[3]] ?: 0) + 1)
                }
            } else {
                firstLine = false
            }
        }

        firstLine = true
        new File(vdjdbSlimPath).splitEachLine("\t") {
            if (!firstLine) {
                if (goodRecord(it) && antigenCountMap[it[3]] >= minCdr3PerAg) {
                    def record
                    recordMap.put(it[1], record = (recordMap[it[0]] ?: new Record(it[0], it[1], it[7], it[9])))
                    record.antigen.add(new AminoAcidSequence(it[3]))
                }
            } else {
                firstLine = false
            }
        }

        recordMap.values()
    }
}
