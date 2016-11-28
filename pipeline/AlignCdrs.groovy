@Grapes(
    @Grab(group='com.milaboratory', module='milib', version='1.1')
)

import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters

import java.util.function.BiFunction

def maxS = args[0], maxID = args[1]

def cdr3AgMap = new HashMap<AminoAcidSequence, Set<String>>()

new File("vdjdb.txt").splitEachLine("\t") { splitLine ->
    cdr3AgMap.compute(new AminoAcidSequence(splitLine[0]), )
}

SequenceTreeMap<AminoAcidSequence, Object> treeMap = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET)

def (mismatches, deletions, insertions, total) = args[1].split(",").collect { it.toInteger() }

def searchParameters = new TreeSearchParameters(mismatches, deletions, insertions, total)

cdr3List.each { treeMap.put(new AminoAcidSequence(it), it) }

new File(args[2]).withPrintWriter { pw -> 
	cdr3List.each {
		def iter = treeMap.getNeighborhoodIterator(new AminoAcidSequence(it), searchParameters)
	
		def to
	
	    while ((to = iter.next()) != null) {
	        if (to != it) {
	            pw.println(iter.currentAlignment.toString().split("[\n ]+")[[1, 4]].join("\t"))
	        }
		}
	}
}
