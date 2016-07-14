package com.antigenomics.cdr3align

import com.antigenomics.vdjdb.scoring.AlignmentScoringProvider
import com.antigenomics.vdjdb.scoring.VdjdbAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import groovyx.gpars.GParsPool
import org.moeaframework.Executor
import org.moeaframework.core.Solution
import org.moeaframework.util.progress.ProgressEvent
import org.moeaframework.util.progress.ProgressListener

import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicInteger

def sout = {
    println "[CDRALING ${new Date().toString()}] $it"
}

def DEFAULT_CONF_THRESHOLD = "1", DEFAULT_SPECIES = "HomoSapiens",
    DEFAULT_GENES = "TRA,TRB", DEFAULT_MIN_CDR3_PER_AG = "2",
    DEFAULT_SCOPE = "5,2,2,7", DEFAULT_PATH = "vdjdb.slim.txt",
    DEFAULT_MOEA_POP_SIZE = "150", DEFAULT_MOEA_GEN = "1000"


def cli = new CliBuilder(usage: "TrainAlignerScoring [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli._(longOpt: "vdjdb-conf-threshold", argName: "0-3", args: 1,
        "VDJdb confidence score threshold. [default=$DEFAULT_CONF_THRESHOLD]")
cli._(longOpt: "species", argName: "HomoSapiens,...", args: 1,
        "Species list. [default=$DEFAULT_SPECIES]")
cli._(longOpt: "genes", argName: "TRA,...", args: 1,
        "Gene list. [default=$DEFAULT_GENES]")
cli._(longOpt: "min-cdr3-per-ag", argName: "1+", args: 1,
        "Minimal number of unique CDR3 sequences per antigen. [default=$DEFAULT_MIN_CDR3_PER_AG]")
cli._(longOpt: "search-scope", argName: "s,i,d,t", args: 1,
        "Overrides CDR3 sequence initial search parameters: allowed number of substitutions (s), insertions (i), " +
                "deletions (d) and total number of mutations. [default=$DEFAULT_SCOPE]")
cli._(longOpt: "vdjdb-slim-path", argName: "path/to/vdjdb.slim.txt", args: 1,
        "Path to vdjdb table in slim format. [default=$DEFAULT_PATH]")
cli._(longOpt: "moea-pop-size", argName: "50-250", args: 1,
        "MOEA population size. [default=$DEFAULT_MOEA_POP_SIZE]")
cli._(longOpt: "moea-gen", argName: "100+", args: 1,
        "Number of MOEA generations to urn. [default=$DEFAULT_MOEA_GEN]")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(2)
}

// requires a pre-built database
// load records

sout "Loading database"

def recordMap = new HashMap<String, Record>()

def antigenCountMap = new HashMap<String, Integer>()

def minDbScore = 2, species = "HomoSapiens"//, gene = "TRB"

def goodRecord = { List<String> record ->
    record[-1].toInteger() >= minDbScore && //record[0] == gene &&
            record[2] == species
}

def firstLine = true
new File("../../database/vdjdb.slim.txt").splitEachLine("\t") {
    if (!firstLine) {
        if (goodRecord(it)) {
            antigenCountMap.put(it[3], (antigenCountMap[it[3]] ?: 0) + 1)
        }
    } else {
        firstLine = false
    }
}

firstLine = true
def minCdr3CountPerAntigen = 5
new File("../../database/vdjdb.slim.txt").splitEachLine("\t") {
    if (!firstLine) {
        if (goodRecord(it) && antigenCountMap[it[3]] >= minCdr3CountPerAntigen) {
            def record
            recordMap.put(it[1], record = (recordMap[it[0]] ?: new Record(it[0], it[1])))
            record.antigen.add(new AminoAcidSequence(it[3]))
        }
    } else {
        firstLine = false
    }
}

sout "Loaded ${recordMap.size()} unique CDR3s."

// align all-vs-all
def searchParameters = new TreeSearchParameters(5, 2, 2, 7)

def treeMap = new SequenceTreeMap<AminoAcidSequence, Record>(AminoAcidSequence.ALPHABET)

recordMap.values().each {
    treeMap.put(it.cdr3, it)
}

sout "Performing alignments."

def alignments = new ConcurrentLinkedQueue<RecordAlignment>()
def counter = new AtomicInteger()

GParsPool.withPool(Runtime.getRuntime().availableProcessors()) {
    treeMap.values().eachParallel { Record from ->
        def iter = treeMap.getNeighborhoodIterator(from.cdr3, searchParameters)
        def to

        def toCdr3Hash = new HashSet<Record>()
        while ((to = iter.next()) != null) {
            if (from.cdr3 != to.cdr3 && from.gene == to.gene && !toCdr3Hash.contains(to)) {
                alignments.add(new RecordAlignment(from, to, iter.currentAlignment))
                toCdr3Hash.add(to)
            }
        }

        int counterVal
        if ((counterVal = counter.incrementAndGet()) % 100 == 0) {
            sout "Done all alignments for $counterVal records, " +
                    "total number of aligned CDR3 pairs is ${alignments.size()}"
        }
    }
}

sout "Done, ${alignments.size()} alignments performed. Proceeding to optimization"

// Run optimization
int popSize = 200, nGenerations = 1000
def listener = new ProgressListener() {
    @Override
    void progressUpdate(ProgressEvent event) {
        sout "[MOEA stats] Number of generations = " + (event.currentNFE / popSize)
    }
}

def problem = new ScoringProblem(alignments)

def result = new Executor()
        .distributeOnAllCores()
        .withProblem(problem)
        .withAlgorithm("NSGAII")
        .withProperty("populationSize", popSize)
        .withMaxEvaluations(nGenerations * popSize)
        .withProgressListener(listener).run()

//display the results
def outputFolder

new File(outputFolder).mkdirs()

def scoringsById = new HashMap<String, VdjdbAlignmentScoring>()
new File("$outputFolder/roc.txt").withPrintWriter { pw ->
    pw.println("id\tsensitivity\tspecificity")
    result.eachWithIndex { Solution solution, int index ->
        pw.println(index + "\t" + solution.getObjective(0) + "\t" + solution.getObjective(1))
        scoringsById.put(index.toString(), problem.decode(solution))
    }
}

AlignmentScoringProvider.saveScoring(scoringsById, "$outputFolder/solutions.txt")