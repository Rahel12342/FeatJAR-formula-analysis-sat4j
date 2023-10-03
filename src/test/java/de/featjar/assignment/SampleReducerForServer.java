package de.featjar.assignment;

import static de.featjar.base.computation.Computations.async;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import de.featjar.base.FeatJAR;
import de.featjar.base.computation.Computations;
import de.featjar.base.computation.IComputation;
import de.featjar.base.io.IO;
import de.featjar.formula.analysis.SampleReducer2;
import de.featjar.formula.analysis.VariableMap;
import de.featjar.formula.analysis.bool.BooleanClauseList;
import de.featjar.formula.analysis.bool.BooleanRepresentationComputation;
import de.featjar.formula.analysis.bool.BooleanSolution;
import de.featjar.formula.analysis.bool.BooleanSolutionList;
import de.featjar.formula.analysis.combinations.LexicographicIterator;
import de.featjar.formula.analysis.sat4j.solver.SAT4JSolutionSolver;
import de.featjar.formula.analysis.sat4j.solver.SAT4JSolver;
import de.featjar.formula.analysis.sat4j.twise.YASA;
import de.featjar.formula.io.FormulaFormats;
import de.featjar.formula.structure.formula.IFormula;
import de.featjar.formula.transformer.ComputeCNFFormula;
import de.featjar.formula.transformer.ComputeNNFFormula;

public class SampleReducerForServer {
	private static VariableMap features = new VariableMap();
	private static BooleanSolutionList fieldVariants = new BooleanSolutionList();
	private static BooleanSolutionList finalFieldVariants = new BooleanSolutionList();
	private static int t;

	public static void main(String[] args) {

		FeatJAR.initialize();
		t = Integer.parseInt(args[0]);

		// initialize feature model & Co for the evaluation tests

		for (int i = 0; i < 5; i++) {
			// Industry1
			IFormula formula = IO.load(Paths.get(
					args[1]),
					FormulaFormats.getInstance()).orElseThrow();

			BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf = Computations.of(formula)
					.map(ComputeNNFFormula::new).map(ComputeCNFFormula::new).map(BooleanRepresentationComputation::new);
			BooleanClauseList clauses = cnf.map(Computations::getKey).computeResult().get();

			fieldVariants = new BooleanSolutionList();
			finalFieldVariants = new BooleanSolutionList();
			features = new VariableMap();

			// get field variants from documents and put them in the finalFieldVariants
			fieldVariants.addAll(getFieldVariants(args[2]));
			for (Iterator<BooleanSolution> iterator = fieldVariants.getAll().iterator(); iterator.hasNext();) {
				BooleanSolution inter = iterator.next();
				int[] array = new int[features.getVariableCount()];
				for (int j = 0; j < features.getVariableCount(); j++) {
					if (inter.size() > j) {
						if (inter.get(j) == 0) {
							array[j] = -(j + 1);
						} else {
							array[j] = j + 1;
						}
					} else {
						array[j] = -(j + 1);
					}
				}
				finalFieldVariants.add(new BooleanSolution(array));
			}
			System.out.println("Number of field configurations: " + finalFieldVariants.size());
			System.out.println("Number of features: " + features.getVariableCount());

			if (i == 0) {
				System.out.println("Sample before reducing");
				checkCoverage(t, Computations.of(clauses), finalFieldVariants, cnf);
				calculateFMInteractions(t, Computations.of(clauses));
				System.out.println();
			}
			List<BooleanSolution> reducedSample = startSampleReducer(clauses, cnf, i, args[3]);
			if (i == 0) {
				System.out.println("Sample after reduction: ");
				checkCoverage(t, Computations.of(clauses), new BooleanSolutionList(reducedSample), cnf);
				System.out.println();
			}
		}
	}

	public static List<BooleanSolution> startSampleReducer(BooleanClauseList clauses,
			BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf, int round, String alg) {
		BooleanSolutionList sampleYASA = new BooleanSolutionList();
		long start = System.currentTimeMillis();

		IComputation<BooleanClauseList> clausesYASA = cnf.map(Computations::getKey);
		sampleYASA = clausesYASA.map(YASA::new).setDependencyComputation(YASA.T, async(t)).compute();
		long end = System.currentTimeMillis();
		float sec = (end - start) / 1000F;
		System.out.println("Seconds needed to calculate the sample of last FM with YASA: " + sec);
		System.out.println("Size of YASA sample for last FM: " + sampleYASA.size());
		System.out.println();

		// reduce the sample
		System.out.println("Reducing Sample: ");
		SampleReducer2 reducer = new SampleReducer2();

		List<BooleanSolution> reducedSample = new ArrayList<>();
		if(alg.equals("R")) {
			reducedSample = reducer.reduceRandom(finalFieldVariants.getAll(), t);
		} else {
			reducedSample = reducer.reduce(finalFieldVariants.getAll(), t);
		}
		System.out.println();

//         reduce the YASA sample to get it ordered by score
		System.out.println("Reducing the YASA sample:");
		List<BooleanSolution> reducedYASASampleList = new ArrayList<>();
		if(alg.equals("R")) {
			reducedSample = reducer.reduceRandom(sampleYASA.getAll(), t);
		} else {
			reducedSample = reducer.reduce(sampleYASA.getAll(), t);
		}
		System.out.println("Sice YASA sample after reduction: " + reducedYASASampleList.size());
		System.out.println();
		if (round == 1) {
			System.out.println("Reverting reduction of reduced sample:");
			reducer.revertReduction(reducedSample, t);
			System.out.println("Reverting reduction of YASA sample:");
			reducer.revertReduction(reducedYASASampleList, t);
		}
		return reducedSample;
	}

	private static ArrayList<BooleanSolution> getFieldVariants(String fieldVariantsFolder) {
		ArrayList<BooleanSolution> clauseList = new ArrayList<>();

		Path root = Paths.get(fieldVariantsFolder);

		// Industry
		File folder = new File(fieldVariantsFolder);
		for (File fileEntry : folder.listFiles()) {

			List<String> lines;
			try {
				lines = Files.readAllLines(root.resolve(fileEntry.getName()));
			} catch (IOException e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}
			lines.remove(0);
			List<Integer> configArray = new ArrayList<>();
			int arrayCounter = 0;
			for (String line : lines) {
				String[] featuresArray = line.split(";");
				if (!featuresArray[0].isEmpty()) {
					if (!configArray.isEmpty()) {
						int[] configArrray = new int[features.getVariableCount() + configArray.size()];
						for (int i = 0; i < configArray.size(); i++) {
							configArrray[i] = configArray.get(i);
						}
						BooleanSolution config = new BooleanSolution(configArrray);
						clauseList.add(config);
						configArray = new ArrayList<>();
						arrayCounter = 0;
					}
				}
				for (int i = 2; i < featuresArray.length; i++) {
					if (featuresArray[i].equals("null")) {
						continue;
					}
					if (features.get(featuresArray[i]).isEmpty()) {
						features.add(featuresArray[i]);
					}
					configArray.add(arrayCounter++, features.get(featuresArray[i]).orElseThrow());
				}
			}
		}

		return clauseList;
	}

	private static void checkCoverage(int t, IComputation<BooleanClauseList> clauses, BooleanSolutionList sample,
			BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf) {
		System.out.println("Sample size:" + sample.size());
		
		// count how many interactions are in the sample overall
		if(t == 1) {
		long countAppearingInteractions = LexicographicIterator.<Void>stream(t, sample.get(0).get().size() * 2)
				.map(interaction -> {
					int[] array = new int[1];
					array[0] = interaction.elementIndices[0] > sample.get(0).get().size() - 1
							? -(interaction.elementIndices[0] - sample.get(0).get().size() + 1)
							: interaction.elementIndices[0] + 1;
					return array;
				}).filter(interaction -> {
					for (BooleanSolution config : sample) {
						if ((interaction[0] < 0 && config.get(-interaction[0] - 1) == interaction[0])
								|| (interaction[0] > 0 && config.get(interaction[0] - 1) == interaction[0])) {
							return true;
						}
					}
					return false;
				}).count();
		System.out.println("Contained interactions in sample:" + countAppearingInteractions);
			
		} else {
			long countAppearingInteractions = LexicographicIterator.<Void>stream(t, sample.get(0).get().size() * 2)
					.map(interaction -> {
						int[] array = new int[2];
						array[0] = interaction.elementIndices[0] > sample.get(0).get().size() - 1
								? -(interaction.elementIndices[0] - sample.get(0).get().size() + 1)
										: interaction.elementIndices[0] + 1;
						array[1] = interaction.elementIndices[1] > sample.get(0).get().size() - 1
								? -(interaction.elementIndices[1] - sample.get(0).get().size() + 1)
										: interaction.elementIndices[1] + 1;
						return array;
					}).filter(interaction -> {
						
						for (BooleanSolution config : sample) {
							if (((interaction[0] < 0 && config.get(-interaction[0] - 1) == interaction[0])
									|| (interaction[0] > 0 && config.get(interaction[0] - 1) == interaction[0]))
									&& ((interaction[1] < 0 && config.get(-interaction[1] - 1) == interaction[1])
											|| (interaction[1] > 0 && config.get(interaction[1] - 1) == interaction[1]))) {
								return true;
							}
						}
						return false;
					}).count();
			System.out.println("Contained interactions in sample:" + countAppearingInteractions);
		}


	}

	public static void calculateFMInteractions(int t, IComputation<BooleanClauseList> clauses) {

		// count how many valid feature interactions exist overall in the feature model
		if(t == 1) {
		long countValidFMInteractions = LexicographicIterator.<Void>stream(t, clauses.compute().getVariableCount() * 2)
				.map(interaction -> {
					int[] array = new int[2];
					array[0] = interaction.elementIndices[0] > clauses.compute().getVariableCount() - 1
							? -(interaction.elementIndices[0] - clauses.compute().getVariableCount() + 1)
							: interaction.elementIndices[0] + 1;
					return array;
				}).filter(interaction -> {
					SAT4JSolver solver = new SAT4JSolutionSolver(clauses.compute());
					if (solver.hasSolution(interaction).get()) {
						return true;
					}

					return false;
				}).count();
		System.out.println("Valid interactions of fm:" + countValidFMInteractions);
		} else {
			long countValidFMInteractions = LexicographicIterator.<Void>stream(t, clauses.compute().getVariableCount() * 2)
					.map(interaction -> {
						int[] array = new int[2];
						array[0] = interaction.elementIndices[0] > clauses.compute().getVariableCount() - 1
								? -(interaction.elementIndices[0] - clauses.compute().getVariableCount() + 1)
										: interaction.elementIndices[0] + 1;
						array[1] = interaction.elementIndices[1] > clauses.compute().getVariableCount() - 1
								? -(interaction.elementIndices[1] - clauses.compute().getVariableCount() + 1)
										: interaction.elementIndices[1] + 1;
						return array;
					}).filter(interaction -> {
						
						SAT4JSolver solver = new SAT4JSolutionSolver(clauses.compute());
						if (solver.hasSolution(interaction).get()) {
							return true;
						}
						
						return false;
					}).count();
			
			System.out.println("Valid interactions of fm:" + countValidFMInteractions);
		}
	}
}
