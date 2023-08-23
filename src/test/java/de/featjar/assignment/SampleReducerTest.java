package de.featjar.assignment;

import static de.featjar.base.computation.Computations.await;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import de.featjar.base.computation.Computations;
import de.featjar.base.computation.IComputation;
import de.featjar.base.io.IO;
import de.featjar.formula.analysis.SampleReducer2;
import de.featjar.formula.analysis.VariableMap;
import de.featjar.formula.analysis.bool.BooleanAssignment;
import de.featjar.formula.analysis.bool.BooleanAssignmentList;
import de.featjar.formula.analysis.bool.BooleanClauseList;
import de.featjar.formula.analysis.bool.BooleanRepresentationComputation;
import de.featjar.formula.analysis.bool.BooleanSolution;
import de.featjar.formula.analysis.bool.BooleanSolutionList;
import de.featjar.formula.analysis.combinations.LexicographicIterator;
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
    	BooleanSolutionList finalFieldVariants = new BooleanSolutionList();

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
        BooleanClauseList clauses = cnf.map(Computations::getKey).computeResult().get();
        
        VariableMap features = Computations.await((IComputation<VariableMap>) cnf.map(Computations::getValue));
        
        fieldVariants.addAll(getFieldVariants(features, "/home/rahel/Dokumente/Repositories/soletta-case-study/020_samples/2015/2015-06-26_18-38-56/random/0"));
        fieldVariants.addAll(getFieldVariants(features, "/home/rahel/Dokumente/Repositories/soletta-case-study/020_samples/2015/2015-06-26_18-38-56/random/0"));

        //        BooleanSolutionList sample = await(clauses.map(YASA::new).set(YASA.T, t));
        
        System.out.println("Base sample");
        
        for(Iterator<BooleanSolution> iterator = fieldVariants.getAll().iterator(); iterator.hasNext();) {
        	BooleanSolution inter = iterator.next();
        	int[] array = new int[inter.size()];
        	for(int i = 0; i < inter.size(); i++) {
        		if(inter.get(i) == 0) {
        			array[i] = -(i+1);
        		} else {
        			array[i] = i+1;
        		}
        	}
        	finalFieldVariants.add(new BooleanSolution(array));
        }
        
        checkCoverage(t, Computations.of(clauses) , finalFieldVariants);

        System.out.println("reduced sample");
        SampleReducer2 reducer = new SampleReducer2();
        List<BooleanSolution> reducedSample = reducer.reduce(null, finalFieldVariants.getAll(), t);
        BooleanSolutionList reducedFieldVariants = new BooleanSolutionList(reducedSample);
        checkCoverage(t, Computations.of(clauses), reducedFieldVariants);

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
		long count = LexicographicIterator.<Void>stream(t, sample.get(0).get().size() * 2).map(interaction -> {int[] array = new int[2];
		array[0] = interaction.elementIndices[0] > sample.get(0).get().size() - 1 ? -(interaction.elementIndices[0] - sample.get(0).get().size() + 1) : interaction.elementIndices[0]+1;
		array[1] = interaction.elementIndices[1] > sample.get(0).get().size() - 1 ? -(interaction.elementIndices[1] - sample.get(0).get().size() + 1) : interaction.elementIndices[1]+1; return array;})
    			.filter(interaction -> {
    		
        	for(BooleanSolution config : sample) {
        		if(((interaction[0] < 0 && config.get(-interaction[0]-1) == interaction[0]) || (interaction[0] > 0 && config.get(interaction[0]-1) == interaction[0])) && 
        				((interaction[1] < 0 && config.get(-interaction[1]-1) == interaction[1]) || (interaction[1] > 0 && config.get(interaction[1]-1) == interaction[1]))) {
        			return true;
        		}
        	}
        	
        	//TODO: valide? Solver aufrufen! dann ist es noch ne extra interaction (um am ende coverage raus zu finden)
        	return false;
    	})
    			.count();	
		System.out.println("Sample size:" + sample.size());
        System.out.println("Sample coverage:" + count);
    }
}
