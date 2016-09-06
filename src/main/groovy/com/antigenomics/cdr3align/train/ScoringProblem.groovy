package com.antigenomics.cdr3align.train

import com.antigenomics.cdr3align.db.RecordAlignment
import com.antigenomics.vdjdb.scoring.VdjdbAlignmentScoring
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import groovy.transform.CompileStatic
import org.moeaframework.core.Solution
import org.moeaframework.core.variable.EncodingUtils
import org.moeaframework.core.variable.RealVariable
import org.moeaframework.problem.AbstractProblem

@CompileStatic
class ScoringProblem extends AbstractProblem {
    final Collection<RecordAlignment> alignments

    static final float MAX_DIAG = 1.0f, MIN_NON_DIAG = -1.0f, MIN_GAP = -1.0f, VAR_FACTOR = 1000f,
                       MAX_SIGMA = 3.0f, SCORE_RANGE = 20f * Math.max(MAX_DIAG, Math.max(-MIN_NON_DIAG, -MIN_GAP))

    static final int N_SUBST = (int) (AminoAcidSequence.ALPHABET.size() * (AminoAcidSequence.ALPHABET.size() + 1) / 2),
                     N_SUBST_1 = AminoAcidSequence.ALPHABET.size(),
                     N_SUBST_2 = N_SUBST_1 * N_SUBST_1,
                     N_VARS = N_SUBST + 1i /*sigma*/ + 1i /*mu*/ + 1i /*gap*/ + 1i /*threshold*/

    ScoringProblem(Collection<RecordAlignment> alignments) {
        super(N_VARS, 2, 1)
        this.alignments = alignments
    }

    @Override
    void evaluate(Solution solution) {
        def scoring = decode(solution)

        int TP = 0, FP = 0, TN = 0, FN = 0, trueExact = 0, totalExact = 0

        alignments.each { RecordAlignment recordAlignment ->
            boolean passThreshold = scoring.computeScore(recordAlignment.alignment) >= scoring.scoreThreshold
            if (recordAlignment.cdr3Match) {
                if (passThreshold) {
                    trueExact++
                }
                totalExact++
            } else {
                if (recordAlignment.antigensMatch) {
                    if (passThreshold) {
                        TP++
                    } else {
                        FN++
                    }
                } else {
                    if (passThreshold) {
                        FP++
                    } else {
                        TN++
                    }
                }
            }
        }

        solution.setConstraint(0, 1.0 - trueExact / (double) totalExact)
        solution.setObjective(0, -TP / (double) Math.max(1, TP + FP))
        solution.setObjective(1, -TP / (double) Math.max(1, TP + FN))
    }

    VdjdbAlignmentScoring decode(Solution solution) {
        def vars = EncodingUtils.getReal(solution)

        def substitutionMatrix = new int[N_SUBST_2]

        int k = 0
        for (int i = 0; i < N_SUBST_1; i++) {
            for (int j = i; j < N_SUBST_1; j++) {
                def var = (int) (VAR_FACTOR * vars[k])
                substitutionMatrix[i * N_SUBST_1 + j] = var
                if (i != j)
                    substitutionMatrix[j * N_SUBST_1 + i] = var
                k++
            }
        }

        def gapPenalty = (int) (VAR_FACTOR * vars[k])
        def scoring = new LinearGapAlignmentScoring(AminoAcidSequence.ALPHABET, substitutionMatrix, gapPenalty)

        def posSigma = (float) vars[++k];
        def posMu = (float) vars[++k];
        def threshold = (float) (VAR_FACTOR * vars[++k])

        new VdjdbAlignmentScoring(scoring, posSigma, posMu, threshold)
    }

    @Override
    Solution newSolution() {
        Solution solution = new Solution(N_VARS, 2, 1)

        int k = 0
        for (int i = 0; i < N_SUBST_1; i++) {
            for (int j = i; j < N_SUBST_1; j++) {
                solution.setVariable(k, i == j ? new RealVariable(0, MAX_DIAG) :
                        new RealVariable(MIN_NON_DIAG, 0))
                k++
            }
        }

        solution.setVariable(k, new RealVariable(MIN_GAP, -2 / VAR_FACTOR))
        solution.setVariable(++k, new RealVariable(0.01, MAX_SIGMA))
        solution.setVariable(++k, new RealVariable(-1.0, 1.0))
        solution.setVariable(++k, new RealVariable(-SCORE_RANGE, SCORE_RANGE))

        solution
    }
}
