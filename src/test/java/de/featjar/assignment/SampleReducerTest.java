package de.featjar.assignment;

import static de.featjar.base.computation.Computations.await;
import static org.junit.jupiter.api.Assertions.assertEquals;

import de.featjar.base.computation.Computations;
import de.featjar.base.computation.IComputation;
import de.featjar.base.io.IO;
import de.featjar.formula.analysis.SampleReducer;
import de.featjar.formula.analysis.bool.BooleanClauseList;
import de.featjar.formula.analysis.bool.BooleanRepresentationComputation;
import de.featjar.formula.analysis.bool.BooleanSolution;
import de.featjar.formula.analysis.bool.BooleanSolutionList;
import de.featjar.formula.analysis.sat4j.twise.CoverageStatistic;
import de.featjar.formula.analysis.sat4j.twise.TWiseStatisticGenerator;
import de.featjar.formula.analysis.sat4j.twise.YASA;
import de.featjar.formula.io.FormulaFormats;
import de.featjar.formula.structure.formula.IFormula;
import de.featjar.formula.transformer.ComputeCNFFormula;
import de.featjar.formula.transformer.ComputeNNFFormula;
import java.nio.file.Paths;
import java.util.List;

public class SampleReducerTest {

    private static final int randomSeed = 1234;

    public static void main(String[] args) {

        IFormula formula = IO.load(
                        Paths.get(
                                "/home/rahel/Dokumente/Repositories/FeatJAR/formula/src/test/resources/GPL/model.xml"),
                        FormulaFormats.getInstance())
                .orElseThrow();
        int t = 2;

        BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf = Computations.of(formula)
                .map(ComputeNNFFormula::new)
                .map(ComputeCNFFormula::new)
                .map(BooleanRepresentationComputation::new);
        IComputation<BooleanClauseList> clauses = cnf.map(Computations::getKey);

        //        BooleanSolutionList sample = await(clauses.map(AllConfigurationGenerator::new));
        BooleanSolutionList sample = await(clauses.map(YASA::new).set(YASA.T, t));
        checkCoverage(t, clauses, sample);

        SampleReducer reducer = new SampleReducer();
        List<BooleanSolution> reducedSample = reducer.reduce(sample.getAll(), t);

        System.out.println();

        //        Logger.logInfo("Reduce...");
        //        long start = System.currentTimeMillis();
        //        List<LiteralList> reducedSample = SampleReducer.reduce(sample, t);
        //        long end = System.currentTimeMillis();
        //
        //        Logger.logInfo("Time: " + ((end - start) / 1000.0));
        //
        //        Collections.sort(solutionList1, Comparator.comparing(LiteralList::toString));
        //        Collections.sort(reducedSample, Comparator.comparing(LiteralList::toString));
        //
        //        System.out.println(solutionList1.size());
        ////        solutionList1.forEach(l -> System.out.println(l.toString()));
        //        System.out.println(reducedSample.size());
        ////        reducedSample.forEach(l -> System.out.println(l.toString()));
    }

    private static void checkCoverage(int t, IComputation<BooleanClauseList> clauses, BooleanSolutionList sample) {
        CoverageStatistic statistic = await(clauses.map(TWiseStatisticGenerator::new)
                .set(TWiseStatisticGenerator.SAMPLE, sample)
                .set(TWiseStatisticGenerator.T, t));
        System.out.println(sample.size());
        assertEquals(1.0, statistic.coverage());
    }
}
