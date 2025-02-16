package com.antigenomics.cdr3align.train

import com.antigenomics.cdr3align.db.Record
import com.antigenomics.cdr3align.db.VdjdbBatchAligner
import com.antigenomics.cdr3align.db.VdjdbLoader
import com.antigenomics.vdjdb.scoring.AlignmentScoringProvider
import com.antigenomics.vdjdb.scoring.VdjdbAlignmentScoring
import org.moeaframework.Executor
import org.moeaframework.core.Solution
import org.moeaframework.util.progress.ProgressEvent
import org.moeaframework.util.progress.ProgressListener

def sout = {
    println "[CDRALING ${new Date().toString()}] $it"
}

def DEFAULT_CONF_THRESHOLD = "1", DEFAULT_SPECIES = "HomoSapiens",
    DEFAULT_GENES = "TRB", DEFAULT_MIN_CDR3_PER_AG = "2",
    DEFAULT_SCOPE = "5,2,2,7", DEFAULT_PATH = "vdjdb.slim.txt", DEFAULT_NEG_PATH = "cdr3align.neg.txt",
    DEFAULT_MOEA_POP_SIZE = "150", DEFAULT_MOEA_GEN = "5000",
    DEFAULT_POS_WEIGHT_LEN = "7",
    DEFAULT_THREADS = Runtime.getRuntime().availableProcessors().toString()

def cli = new CliBuilder(usage: "TrainAlignerScoring [options] output_folder/")
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
cli._(longOpt: "pos-weight-len", argName: "1..10", args: 1,
        "Length of positional weight vector (should be odd number). [default=$DEFAULT_POS_WEIGHT_LEN]")
cli._(longOpt: "vdjdb-slim-path", argName: "path/to/vdjdb.slim.txt", args: 1,
        "Path to vdjdb table in slim format. [default=$DEFAULT_PATH]")
cli._(longOpt: "neg-control-path", argName: "cdr3align.neg.txt", args: 1,
        "Path to table with negative control set (columns gene | cdr3 | v | j). [default=$DEFAULT_NEG_PATH]")
cli._(longOpt: "only-v-match", "Only evaluate alignments between CDR3 sequences " +
        "that have at least one matching V segment.")
cli._(longOpt: "only-j-match", "Only evaluate alignments between CDR3 sequences " +
        "that have at least one matching V segment.")
cli._(longOpt: "moea-pop-size", argName: "50-250", args: 1,
        "MOEA population size. [default=$DEFAULT_MOEA_POP_SIZE]")
cli._(longOpt: "moea-gen", argName: "100+", args: 1,
        "Number of MOEA generations to run. [default=$DEFAULT_MOEA_GEN]")
cli._(longOpt: "threads", argName: "1+", args: 1,
        "Number of threads to distribute on. [default=$DEFAULT_THREADS]")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(1)
}

if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(1)
}

def outputFolder = opt.arguments()[-1],
    vdjdbConfThreshold = (opt.'vdjdb-conf-threshold' ?: DEFAULT_CONF_THRESHOLD).toInteger(),
    species = (opt.'species' ?: DEFAULT_SPECIES).split(",") as List<String>,
    genes = (opt.'genes' ?: DEFAULT_GENES).split(",") as List<String>,
    minCdr3PerAg = (opt.'min-cdr3-per-ag' ?: DEFAULT_MIN_CDR3_PER_AG).toInteger(),
    searchScope = (opt.'search-scope' ?: DEFAULT_SCOPE).split(",").collect { it.toInteger() },
    posWeightLen = (opt.'pos-weight-len' ?: DEFAULT_POS_WEIGHT_LEN).toInteger(),
    vdjdbSlimPath = opt.'vdjdb-slim-path' ?: DEFAULT_PATH,
    negControlPath = opt.'neg-control-path' ?: DEFAULT_NEG_PATH,
    onlyVMatch = (boolean) opt.'only-v-match', onlyJMatch = (boolean) opt.'only-j-match',
    moeaPopSize = (opt.'moea-pop-size' ?: DEFAULT_MOEA_POP_SIZE).toInteger(),
    moeaGen = (opt.'moea-gen' ?: DEFAULT_MOEA_GEN).toInteger(),
    threads = (opt.'threads' ?: DEFAULT_THREADS).toInteger()

// requires a pre-built database
// load records

sout "Loading database.."

def records = new VdjdbLoader(vdjdbSlimPath).load(vdjdbConfThreshold, genes, species, minCdr3PerAg) as List

sout "Loaded ${records.size()} unique CDR3s."

if (new File(negControlPath).exists()) {
    sout "Loading negative set.."

    def existingCdr3 = new HashSet<String>()
    records.each { existingCdr3.add(it.cdr3.toString()) }

    new File(negControlPath).splitEachLine("\t") { splitLine ->
        def gene = splitLine[0], cdr3 = splitLine[1]

        if (genes.any { it.equalsIgnoreCase(gene) } && !existingCdr3.contains(cdr3)) { // prevent exact matches
            records.add(new Record(gene, cdr3, splitLine[2], splitLine[3]))
        }
    }

    sout "Loaded in total ${records.size()} unique CDR3s."
}

// align all-vs-all
def batchAligner = new VdjdbBatchAligner(searchScope[0], searchScope[2], searchScope[1], searchScope[3],
        onlyVMatch, onlyJMatch)

batchAligner.add(records)

sout "Performing alignments."

def alignments = batchAligner.align(threads)

sout "Done, ${alignments.size()} alignments performed. Proceeding to optimization"

// Run optimization
def listener = new ProgressListener() {
    @Override
    void progressUpdate(ProgressEvent event) {
        int ngen = event.currentNFE / moeaPopSize
        if (ngen % 10 == 0) {
            sout "[MOEA stats] Number of generations = " + ngen
        }
    }
}

def problem = new ScoringProblem(alignments, posWeightLen)

def result = new Executor()
        .distributeOn(threads)
        .withProblem(problem)
        .withAlgorithm("NSGAII")
        .withProperty("populationSize", moeaPopSize)
        .withMaxEvaluations(moeaGen * moeaPopSize)
        .withProgressListener(listener).run()

sout "Finished. Writing results."

new File(outputFolder).mkdirs()

def scoringsById = new HashMap<String, VdjdbAlignmentScoring>()
new File("$outputFolder/roc.txt").withPrintWriter { pw ->
    pw.println("id\tprecision\trecall\texact.misclassified")
    result.eachWithIndex { Solution solution, int index ->
        pw.println(index + "\t" +
                (-solution.getObjective(0)) + "\t" +
                (-solution.getObjective(1)) + "\t" +
                (solution.getConstraint(0)))
        scoringsById.put(index.toString(), problem.decode(solution))
    }
}

AlignmentScoringProvider.saveScoring(scoringsById, "$outputFolder/solutions.txt")