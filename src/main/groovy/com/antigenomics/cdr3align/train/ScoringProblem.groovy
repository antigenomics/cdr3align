package com.antigenomics.cdr3align.train

import com.antigenomics.cdr3align.db.RecordAlignment
import com.antigenomics.vdjdb.scoring.VdjdbAlignmentScoring
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import org.moeaframework.core.Solution
import org.moeaframework.core.variable.EncodingUtils
import org.moeaframework.core.variable.RealVariable
import org.moeaframework.problem.AbstractProblem

class ScoringProblem extends AbstractProblem {
    final Collection<RecordAlignment> alignments
    final int nPositionalWeights

    static final float MAX_DIAG = 1.0, MIN_NON_DIAG = -1.0, MIN_GAP = -1.0, VAR_FACTOR = 1000,
                       SCORE_RANGE = 20 * Math.max(MAX_DIAG, Math.max(-MIN_NON_DIAG, -MIN_GAP))

    static final int N_SUBST = AminoAcidSequence.ALPHABET.size() * (AminoAcidSequence.ALPHABET.size() + 1) / 2,
                     N_SUBST_1 = AminoAcidSequence.ALPHABET.size(),
                     N_SUBST_2 = N_SUBST_1 * N_SUBST_1,
                     N_VARS = N_SUBST + 1 /*gap*/ + 1 /*threshold*/,
                     N_OBJ = 3

    ScoringProblem(Collection<RecordAlignment> alignments,
                   int nPositionalWeights = 7) {
        super(N_VARS + nPositionalWeights, N_OBJ)
        this.alignments = alignments
        this.nPositionalWeights = nPositionalWeights
        assert nPositionalWeights % 2 == 1
    }

    @Override
    void evaluate(Solution solution) {
        def scoring = decode(solution)

        int TP = 0, FP = 0, TN = 0, FN = 0, trueExact = 0, totalExact = 0

        alignments.each { RecordAlignment recordAlignment ->
            double score = scoring.computeScore(recordAlignment.alignment)
            if (recordAlignment.tcrData1 == recordAlignment.tcrData2) {
                if (score >= scoring.scoreThreshold) {
                    trueExact++
                }
                totalExact++
            } else {
                if (recordAlignment.antigensMatch) {
                    if (score >= scoring.scoreThreshold) {
                        TP++
                    } else {
                        FN++
                    }
                } else {
                    if (score >= scoring.scoreThreshold) {
                        FP++
                    } else {
                        TN++
                    }
                }
            }
        }

        solution.setObjective(0, -TP / (double) Math.max(1, TP + FP))
        solution.setObjective(1, -TP / (double) Math.max(1, TP + FN))
        solution.setObjective(2, -trueExact / (double) totalExact)
    }

    VdjdbAlignmentScoring decode(Solution solution) {
        double[] vars = EncodingUtils.getReal(solution)

        int[] substitutionMatrix = new int[N_SUBST_2]

        int k = 0
        for (int i = 0; i < N_SUBST_1; i++) {
            for (int j = i; j < N_SUBST_1; j++) {
                int var = VAR_FACTOR * vars[k]
                substitutionMatrix[i * N_SUBST_1 + j] = var
                if (i != j)
                    substitutionMatrix[j * N_SUBST_1 + i] = var
                k++
            }
        }

        int gapPenalty = VAR_FACTOR * vars[k]
        def scoring = new LinearGapAlignmentScoring(AminoAcidSequence.ALPHABET, substitutionMatrix, gapPenalty)

        float[] posWeights = new float[nPositionalWeights]
        for (int i = 0; i < nPositionalWeights; i++) {
            posWeights[i] = vars[++k]
        }
        def threshold = (float) (VAR_FACTOR * vars[++k])

        new VdjdbAlignmentScoring(scoring, posWeights, threshold)
    }

    @Override
    Solution newSolution() {
        Solution solution = new Solution(N_VARS + nPositionalWeights, N_OBJ)

        int k = 0
        for (int i = 0; i < N_SUBST_1; i++) {
            for (int j = i; j < N_SUBST_1; j++) {
                solution.setVariable(k, i == j ? new RealVariable(0, MAX_DIAG) :
                        new RealVariable(MIN_NON_DIAG, 0))
                k++
            }
        }

        solution.setVariable(k, new RealVariable(MIN_GAP, -2 / VAR_FACTOR))

        for (int i = 0; i < nPositionalWeights; i++) {
            solution.setVariable(++k, new RealVariable(0, 1))
        }

        solution.setVariable(++k, new RealVariable(-SCORE_RANGE, SCORE_RANGE))

        solution
    }
}
