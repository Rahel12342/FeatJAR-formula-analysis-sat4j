package de.featjar.assignment;

import static de.featjar.base.computation.Computations.await;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import de.featjar.base.computation.Computations;
import de.featjar.base.computation.IComputation;
import de.featjar.base.io.IO;
import de.featjar.formula.analysis.SampleReducer2;
import de.featjar.formula.analysis.VariableMap;
import de.featjar.formula.analysis.bool.BooleanClauseList;
import de.featjar.formula.analysis.bool.BooleanRepresentationComputation;
import de.featjar.formula.analysis.bool.BooleanSolution;
import de.featjar.formula.analysis.bool.BooleanSolutionList;
import de.featjar.formula.analysis.sat4j.twise.CoverageStatistic;
import de.featjar.formula.analysis.sat4j.twise.TWiseStatisticGenerator;
import de.featjar.formula.io.FormulaFormats;
import de.featjar.formula.structure.formula.IFormula;
import de.featjar.formula.transformer.ComputeCNFFormula;
import de.featjar.formula.transformer.ComputeNNFFormula;

public class SampleReducerTest {

    private static final int randomSeed = 1234;

    public static void main(String[] args) {
    	BooleanSolutionList fieldVariants = new BooleanSolutionList();

        IFormula formula = IO.load(
                        Paths.get(
                                "/home/rahel/Dokumente/Repositories/soletta-case-study/010_models/2015/2015-06-26_18-38-56/model.xml"),
                        FormulaFormats.getInstance())
                .orElseThrow();
        int t = 2;
        
        
        
        BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf = Computations.of(formula)
                .map(ComputeNNFFormula::new)
                .map(ComputeCNFFormula::new)
                .map(BooleanRepresentationComputation::new);
        IComputation<BooleanClauseList> clauses = cnf.map(Computations::getKey);
        
        VariableMap features = Computations.await((IComputation<VariableMap>) cnf.map(Computations::getValue));
        
        fieldVariants.addAll(getFieldVariants(features, "/home/rahel/Dokumente/Repositories/soletta-case-study/020_samples/2015/2015-06-26_18-38-56/random/0"));
        fieldVariants.addAll(getFieldVariants(features, "/home/rahel/Dokumente/Repositories/soletta-case-study/020_samples/2015/2015-06-26_18-38-56/random/0"));
//
//        //        BooleanSolutionList sample = await(clauses.map(AllConfigurationGenerator::new));
//        BooleanSolutionList sample = await(clauses.map(YASA::new).set(YASA.T, t));

        checkCoverage(t, clauses, fieldVariants);

        SampleReducer2 reducer = new SampleReducer2();
        List<BooleanSolution> reducedSample = reducer.reduce(null, fieldVariants.getAll(), t);
        BooleanSolutionList reducedFieldVariants = new BooleanSolutionList(reducedSample);
        checkCoverage(t, clauses, reducedFieldVariants);

//        System.out.println(reducedSample.size());

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

	private static ArrayList<BooleanSolution> getFieldVariants(VariableMap features, String fieldVariantsFolder) {

		File folder = new File(fieldVariantsFolder);
		
		ArrayList<BooleanSolution> clauseList = new ArrayList<>();
		
		Path root = Paths.get(fieldVariantsFolder);
		for (File fileEntry : folder.listFiles()) {

			List<String> lines;
			try {
				lines = Files.readAllLines(root.resolve(fileEntry.getName()));
			} catch (IOException e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}
			int[] configArray = new int[features.getVariableCount()];
			int arrayCounter = 0;
			for (String line : lines) {
				configArray[arrayCounter++] = features.get(line).orElseThrow();
			}


			BooleanSolution config = new BooleanSolution(configArray);
			clauseList.add(config);
		}
		return clauseList;
	}

	private static void checkCoverage(int t, IComputation<BooleanClauseList> clauses, BooleanSolutionList sample) {
        CoverageStatistic statistic = await(clauses.map(TWiseStatisticGenerator::new)
                .set(TWiseStatisticGenerator.SAMPLE, sample)
                .set(TWiseStatisticGenerator.T, t));
        System.out.println(sample.size());
        System.out.println(statistic.coverage());
//        assertEquals(1.0, statistic.coverage());
    }
}
