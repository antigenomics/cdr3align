import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.mutations.MutationType;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;
import com.milaboratory.core.tree.TreeSearchParameters;

public class AlignCdrAux {
    public static void main(String[] args) throws IOException {
        int maxSubstitutions = Integer.parseInt(args[0]),
                maxIndels = Integer.parseInt(args[1]);

        String inputFileName = args[2], outputFileName = args[3];

        Map<AminoAcidSequence, Cdr3Info> cdr3AntigenMap = new HashMap<>();
        try (Stream<String> stream = Files.lines(new File(inputFileName).toPath())) {
            stream.forEach(line -> {
                String[] splitString = line.split("\t");
                cdr3AntigenMap.compute(new AminoAcidSequence(splitString[0]),
                        (aminoAcidSequence, cdr3Info) -> {
                            if (cdr3Info == null) {
                                cdr3Info = new Cdr3Info(aminoAcidSequence);
                            }
                            cdr3Info.addAntigen(splitString[1]);
                            return cdr3Info;
                        });
            });
        }

        System.out.println("Loaded " + cdr3AntigenMap.size() + " cdr3 sequences.");

        SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);

        cdr3AntigenMap.entrySet().forEach(kvp -> stm.put(kvp.getKey(), kvp.getValue()));

        String header = "substs\tindels\t" +
                "same.ag\tcdr3.len\tweight\t" +
                "mut.type\tmut.pos\tmut.from\tmut.to";


        PrintWriter pw = new PrintWriter(new File(outputFileName));

        pw.println(header);

        for (int i = 0; i <= maxSubstitutions; i++) {
            for (int j = 0; j <= maxIndels; j++) {
                TreeSearchParameters tsp = new TreeSearchParameters(i, j, j);
                String paramString = i + "\t" + j + "\t";

                Queue<String> lines = new ConcurrentLinkedQueue<>();

                AtomicInteger counter = new AtomicInteger(), mutationCounter = new AtomicInteger();

                int ii = i;
                int jj = j;
                cdr3AntigenMap.values().parallelStream().forEach(thisCdr3Info -> {
                            Cdr3Info otherCdr3Info;
                            NeighborhoodIterator iter = stm.getNeighborhoodIterator(thisCdr3Info.cdr3, tsp);

                            Map<AminoAcidSequence, List<AlignmentInfo>> alignmentVariants = new HashMap<>();

                            while ((otherCdr3Info = (Cdr3Info) iter.next()) != null) {
                                if (thisCdr3Info.nonDuplicateComparison(otherCdr3Info)) {
                                    AlignmentInfo alignmentInfo = new AlignmentInfo(
                                            thisCdr3Info.antigensOverlap(otherCdr3Info),
                                            iter.getCurrentAlignment());

                                    alignmentVariants.compute(otherCdr3Info.cdr3,
                                            (cdr3, alignments) -> {
                                                if (alignments == null) {
                                                    alignments = new ArrayList<>();
                                                }
                                                alignments.add(alignmentInfo);
                                                return alignments;
                                            });
                                }
                            }

                            int mCount = 0;

                            for (List<AlignmentInfo> alignmentInfos : alignmentVariants.values()) {
                                float weight = 1.0f / alignmentInfos.size();

                                for (AlignmentInfo alignmentInfo : alignmentInfos) {
                                    Alignment alignment = alignmentInfo.alignment;

                                    String prefix = paramString + (alignmentInfo.sameAntigen ? 1 : 0) + "\t" +
                                            alignment.getSequence1().size() + "\t" + weight + "\t";

                                    Mutations mutations = alignment.getAbsoluteMutations();

                                    for (int k = 0; k < mutations.size(); k++) {
                                        MutationType mutationType = mutations.getTypeByIndex(k);

                                        lines.add(prefix +
                                                shortMutationType(mutationType) + "\t" +
                                                mutations.getPositionByIndex(k) + "\t" +
                                                (mutationType == MutationType.Insertion ?
                                                        mutations.getFromAsSymbolByIndex(k) : "-") + "\t" +
                                                (mutationType == MutationType.Deletion ?
                                                        mutations.getToAsSymbolByIndex(k) : "-")
                                        );

                                        mCount = mutationCounter.incrementAndGet();
                                    }
                                }
                            }

                            int count = counter.incrementAndGet();

                            if (count % 100 == 0) {
                                System.out.println("[" + (new Date()) + "] substs=" + ii + " indels=" + jj +
                                        " queried " + count + " of " + cdr3AntigenMap.size() + " cdr3 sequences. " +
                                        "Recorded " + mCount + " mutations so far.");
                            }
                        }
                );

                lines.forEach(pw::println);
            }
        }

        pw.close();
    }

    private static class Cdr3Info {
        final AminoAcidSequence cdr3;
        final Set<String> antigens = new HashSet<>();

        Cdr3Info(AminoAcidSequence cdr3) {
            this.cdr3 = cdr3;
        }

        void addAntigen(String antigen) {
            antigens.add(antigen);
        }

        boolean antigensOverlap(Cdr3Info other) {
            return overlaps(this.antigens, other.antigens);
        }

        boolean nonDuplicateComparison(Cdr3Info otherCdr3Info) {
            return cdr3.compareTo(otherCdr3Info.cdr3) > 0;
        }

        private static boolean overlaps(Set<String> set1, Set<String> set2) {
            if (set1.size() > set2.size()) {
                return overlaps(set2, set1);
            }
            for (String value : set1) {
                if (set2.contains(value))
                    return true;
            }
            return false;
        }
    }

    private static class AlignmentInfo {
        final boolean sameAntigen;
        final Alignment alignment;

        AlignmentInfo(boolean sameAntigen, Alignment alignment) {
            this.sameAntigen = sameAntigen;
            this.alignment = alignment;
        }
    }

    static String shortMutationType(MutationType mutationType) {
        switch (mutationType) {
            case Substitution:
                return "S";
            case Insertion:
                return "I";
            case Deletion:
                return "D";
            default:
                throw new IllegalArgumentException();
        }
    }
}
