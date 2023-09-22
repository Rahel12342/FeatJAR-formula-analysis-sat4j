package de.featjar.assignment;

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
import de.featjar.formula.io.FormulaFormats;
import de.featjar.formula.structure.formula.IFormula;
import de.featjar.formula.transformer.ComputeCNFFormula;
import de.featjar.formula.transformer.ComputeNNFFormula;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class SampleReducerTest {

    private static VariableMap features = new VariableMap();
    private static BooleanSolutionList fieldVariants = new BooleanSolutionList();
    private static BooleanSolutionList finalFieldVariants = new BooleanSolutionList();

    public static void main(String[] args) {

        FeatJAR.initialize();
        int t = 1;
        // initialize feature model & Co for the evaluation tests

        for (int i = 0; i < 5; i++) {
            //			//Industry1
            //			IFormula formula = IO.load(Paths.get(
            //					"/home/rahel/Dokumente/Daten MA
            // Industry/fieldsamplingdata/data/industry/1/2023-07-24/obfuscated.xml"),
            //					FormulaFormats.getInstance()).orElseThrow();
            //
            //			//Industry2
            //			IFormula formula = IO.load(Paths.get(
            //					"/home/rahel/Dokumente/Daten MA
            // Industry/fieldsamplingdata/data/industry/2/2023-07-24/obfuscated.xml"),
            //					FormulaFormats.getInstance()).orElseThrow();
            //
            //			//Industry3
            //			IFormula formula = IO.load(Paths.get(
            //					"/home/rahel/Dokumente/Daten MA
            // Industry/fieldsamplingdata/data/industry/3/2023-07-24/obfuscated.xml"),
            //					FormulaFormats.getInstance()).orElseThrow();
            //
            //			//Industry4
            //			IFormula formula = IO.load(Paths.get(
            //					"/home/rahel/Dokumente/Daten MA
            // Industry/fieldsamplingdata/data/industry/4/2023-07-24/obfuscated.xml"),
            //					FormulaFormats.getInstance()).orElseThrow();
            //
            //			//Industry5
            //			IFormula formula = IO.load(Paths.get(
            //					"/home/rahel/Dokumente/Daten MA
            // Industry/fieldsamplingdata/data/industry/5/2023-07-24/obfuscated.xml"),
            //					FormulaFormats.getInstance()).orElseThrow();

            //			 Soletta Case Study
            //			IFormula formula = IO.load(Paths.get(
            //
            //	"/home/rahel/Dokumente/Repositories/soletta-case-study/010_models/2017/2017-03-09_21-02-40/model.xml"),
            //					FormulaFormats.getInstance()).orElseThrow();

            //			 Soletta Case Study
            //			IFormula formula = IO.load(Paths.get(
            //
            //	"/home/rahel/Dokumente/Repositories/soletta-case-study/010_models/2017/2017-03-09_21-02-40/model.xml"),
            //					FormulaFormats.getInstance()).orElseThrow();

            //			 ERP
            IFormula formula = IO.load(
                            Paths.get(
                                    "/home/rahel/Dokumente/Repositories/guided-config-evo-eval-data/quantitative_evaluation/from_product_engineers_perspective/evaluation_software_and_data/ERP/ERP.xml"),
                            FormulaFormats.getInstance())
                    .orElseThrow();

            BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf = Computations.of(formula)
                    .map(ComputeNNFFormula::new)
                    .map(ComputeCNFFormula::new)
                    .map(BooleanRepresentationComputation::new);
            BooleanClauseList clauses =
                    cnf.map(Computations::getKey).computeResult().get();

            fieldVariants = new BooleanSolutionList();
            finalFieldVariants = new BooleanSolutionList();
            //			features = Computations.await((IComputation<VariableMap>) cnf.map(Computations::getValue));
            features = new VariableMap();

            // get field variants from documents and put them in the finalFieldVariants
            fieldVariants.addAll(getFieldVariantsComplete());
            for (Iterator<BooleanSolution> iterator = fieldVariants.getAll().iterator(); iterator.hasNext(); ) {
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

            System.out.println("Number of features: " + features.getVariableCount());

            //			BooleanSolutionList sampleYASA = new BooleanSolutionList();
            //			IComputation<BooleanClauseList> clausesYASA = cnf.map(Computations::getKey);
            //			sampleYASA = clausesYASA.map(YASA::new).setDependencyComputation(YASA.T, async(t)).compute();
            //			for (Iterator<BooleanSolution> iterator = sampleYASA.getAll().iterator(); iterator.hasNext();) {
            //				BooleanSolution inter = iterator.next();
            //				int[] array = new int[features.getVariableCount()];
            //				for (int j = 0; j < features.getVariableCount(); j++) {
            //					if (inter.size() > j) {
            //						if (inter.get(j) == 0) {
            //							array[j] = -(j + 1);
            //						} else {
            //							array[j] = j + 1;
            //						}
            //					} else {
            //						array[j] = -(j + 1);
            //					}
            //				}
            //				finalFieldVariants.add(new BooleanSolution(array));
            //			}
            //			if(i == 0) {
            //				System.out.println("Sample before reducing");
            //				checkCoverage(t, Computations.of(clauses), finalFieldVariants, cnf);
            //				calculateFMInteractions(t, Computations.of(clauses));
            ////				checkCoverage(t, null, finalFieldVariants, null);
            //				System.out.println();
            //			}
            List<BooleanSolution> reducedSample = startSampleReducer(clauses, cnf, i);
            //			List<BooleanSolution> reducedSample = startSampleReducer(null, null, i);
            //			if (i == 0) {
            //				System.out.println("Sample after reduction: ");
            //				checkCoverage(t, Computations.of(clauses), new BooleanSolutionList(reducedSample), cnf);
            ////				checkCoverage(t, null, new BooleanSolutionList(reducedSample), null);
            //				System.out.println();
            //			}
        }
    }

    public static List<BooleanSolution> startSampleReducer(
            BooleanClauseList clauses, BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf, int round) {
        int t = 1;

        BooleanSolutionList sampleYASA = new BooleanSolutionList();
        long start = System.currentTimeMillis();

        //		 calculate sample of last fm with YASA
        //		IComputation<BooleanClauseList> clausesYASA = cnf.map(Computations::getKey);
        //		sampleYASA = clausesYASA.map(YASA::new).setDependencyComputation(YASA.T, async(t)).compute();
        //		long end = System.currentTimeMillis();
        //		float sec = (end - start) / 1000F;
        //		System.out.println("Seconds needed to calculate the sample of last FM with YASA: " + sec);
        //		System.out.println("Size of YASA sample for last FM: " + sampleYASA.size());
        //		System.out.println();

        // reduce the sample
        System.out.println("Reducing Sample: ");
        SampleReducer2 reducer = new SampleReducer2();

        List<BooleanSolution> reducedSample = reducer.reduceRandom(finalFieldVariants.getAll(), t);
        System.out.println();

        // reduce the YASA sample to get it ordered by score
        //		System.out.println("Reducing the YASA sample:");
        //		List<BooleanSolution> reducedYASASampleList = reducer.reduce(sampleYASA.getAll(), t);
        //		System.out.println("Sice YASA sample after reduction: " + reducedYASASampleList.size());
        //		System.out.println();
        //		if(round == 0) {
        //			System.out.println("Reverting reduction of reduced sample:");
        //			reducer.revertReduction(reducedSample, t);
        ////			System.out.println("Reverting reduction of YASA sample:");
        ////			reducer.revertReduction(reducedYASASampleList, t);
        //		}
        return reducedSample;
    }

    private static BooleanSolutionList getFieldVariantsComplete() {
        BooleanSolutionList fieldVariants = new BooleanSolutionList();

        // Agrib
        //		fieldVariants.addAll(getFieldVariants(
        //
        //	"/home/rahel/Dokumente/Repositories/guided-config-evo-eval-data/quantitative_evaluation/from_product_engineers_perspective/evaluation_software_and_data/Agri"));

        // ERP
        fieldVariants.addAll(
                getFieldVariants(
                        "/home/rahel/Dokumente/Repositories/guided-config-evo-eval-data/quantitative_evaluation/from_product_engineers_perspective/evaluation_software_and_data/ERP"));

        // industry partner

        // soletta case study
        //		String startPath = "/home/rahel/Dokumente/Repositories/soletta-case-study/020_samples/2015/";
        //		String endPath = "/random/0";
        //
        //		File folder = new File(startPath);
        //
        //		for (File fileEntry : folder.listFiles()) {
        //			String file = fileEntry.getName();
        //			fieldVariants.addAll(getFieldVariants(startPath + file + endPath));
        //		}
        //
        //		startPath = "/home/rahel/Dokumente/Repositories/soletta-case-study/020_samples/2016/";
        //		endPath = "/random/0";
        //
        //		folder = new File(startPath);
        //
        //		for (File fileEntry : folder.listFiles()) {
        //			String file = fileEntry.getName();
        //			fieldVariants.addAll(getFieldVariants(startPath + file + endPath));
        //		}
        //
        //		startPath = "/home/rahel/Dokumente/Repositories/soletta-case-study/020_samples/2017/";
        //		endPath = "/random/0";
        //
        //		folder = new File(startPath);
        //
        //		for (File fileEntry : folder.listFiles()) {
        //			String file = fileEntry.getName();
        //			fieldVariants.addAll(getFieldVariants(startPath + file + endPath));
        //		}

        return fieldVariants;
    }

    private static ArrayList<BooleanSolution> getFieldVariants(String fieldVariantsFolder) {
        ArrayList<BooleanSolution> clauseList = new ArrayList<>();

        Path root = Paths.get(fieldVariantsFolder);

        // Agrib & ERP
        File folder = new File(fieldVariantsFolder + "/Config.csv");
        List<String> lines = new ArrayList<>();
        try {
            lines = Files.readAllLines(root.resolve(folder.getName()));
        } catch (IOException e) {
            e.printStackTrace();
        }

        int configID = 0;
        List<String> configuration = new ArrayList<>();

        for (int i = 0; i < lines.size(); i++) {
            String[] entry = lines.get(i).split(",");
            int thisConfigID = Integer.parseInt(entry[0]);
            String thisFeature = entry[1];
            if (i == 0) {
                configID = thisConfigID;
            }
            if (features.get(thisFeature).isEmpty()) {
                features.add(thisFeature);
            }
            if (configID == thisConfigID) {
                configuration.add(thisFeature);
            } else {
                configID = thisConfigID;
                int[] configArray = new int[configuration.size() + features.getVariableCount()];
                for (int j = 0; j < configuration.size(); j++) {
                    configArray[j] = features.get(thisFeature).get();
                }
                BooleanSolution config = new BooleanSolution(configArray);
                clauseList.add(config);
                configuration = new ArrayList<>();
                configuration.add(thisFeature);
            }
        }

        //		//Soletta Case Study
        //		File folder = new File(fieldVariantsFolder);
        //		for (File fileEntry : folder.listFiles()) {
        //
        //			List<String> lines;
        //			try {
        //				lines = Files.readAllLines(root.resolve(fileEntry.getName()));
        //			} catch (IOException e) {
        //				e.printStackTrace();
        //				throw new RuntimeException(e);
        //			}
        //			int[] configArray = new int[lines.size() + features.getVariableCount()];
        //			int arrayCounter = 0;
        //			for (String line : lines) {
        //				if (features.get(line).isEmpty()) {
        //					features.add(line);
        //				}
        //				configArray[arrayCounter++] = features.get(line).orElseThrow();
        //			}
        //
        //			BooleanSolution config = new BooleanSolution(configArray);
        //			clauseList.add(config);
        //		}

        return clauseList;
    }

    private static void checkCoverage(
            int t,
            IComputation<BooleanClauseList> clauses,
            BooleanSolutionList sample,
            BooleanRepresentationComputation<IFormula, BooleanClauseList> cnf) {
        System.out.println("Sample size:" + sample.size());

        // t = 1
        // count how many interactions are in the sample overall
        long countAppearingInteractions = LexicographicIterator.<Void>stream(
                        t, sample.get(0).get().size() * 2)
                .map(interaction -> {
                    int[] array = new int[1];
                    array[0] =
                            interaction.elementIndices[0] > sample.get(0).get().size() - 1
                                    ? -(interaction.elementIndices[0]
                                            - sample.get(0).get().size()
                                            + 1)
                                    : interaction.elementIndices[0] + 1;
                    return array;
                })
                .filter(interaction -> {
                    for (BooleanSolution config : sample) {
                        if ((interaction[0] < 0 && config.get(-interaction[0] - 1) == interaction[0])
                                || (interaction[0] > 0 && config.get(interaction[0] - 1) == interaction[0])) {
                            return true;
                        }
                    }
                    return false;
                })
                .count();
        System.out.println("Contained interactions in sample:" + countAppearingInteractions);

        // count how many of the interactions in the sample also appear in the feature
        // model
        //		long countValidToFmSampleInteractions = LexicographicIterator.<Void>stream(t, sample.get(0).get().size() *
        // 2)
        //				.map(interaction -> {
        //					int[] array = new int[2];
        //					array[0] = interaction.elementIndices[0] > sample.get(0).get().size() - 1
        //							? -(interaction.elementIndices[0] - sample.get(0).get().size() + 1)
        //							: interaction.elementIndices[0] + 1;
        //					return array;
        //				}).filter(interaction -> {
        //
        //					SAT4JSolver solver = new SAT4JSolutionSolver(clauses.compute());
        //					VariableMap fMfeatures = Computations
        //							.await((IComputation<VariableMap>) cnf.map(Computations::getValue));
        //
        //					Result<String> fmFeatureName1 = features.get(Math.abs(interaction[0]));
        //
        //					if (!fMfeatures.get(fmFeatureName1.get()).isEmpty()) {
        //						for (BooleanSolution config : sample) {
        //							if ((interaction[0] < 0 && config.get(-interaction[0] - 1) == interaction[0])
        //									|| (interaction[0] > 0 && config.get(interaction[0] - 1) == interaction[0])) {
        //								int[] featureInteractions = new int[1];
        //								featureInteractions[0] = interaction[0] > 0 ? fMfeatures.get(fmFeatureName1.get()).get()
        //										: -fMfeatures.get(fmFeatureName1.get()).get();
        //
        //								if (solver.hasSolution(featureInteractions).get()) {
        //									return true;
        //								}
        //							}
        //						}
        //					}
        //
        //					return false;
        //				}).count();

        // t = 2
        // count how many interactions are in the sample overall
        //		long countAppearingInteractions = LexicographicIterator.<Void>stream(t, sample.get(0).get().size() * 2)
        //				.map(interaction -> {
        //					int[] array = new int[2];
        //					array[0] = interaction.elementIndices[0] > sample.get(0).get().size() - 1
        //							? -(interaction.elementIndices[0] - sample.get(0).get().size() + 1)
        //							: interaction.elementIndices[0] + 1;
        //					array[1] = interaction.elementIndices[1] > sample.get(0).get().size() - 1
        //							? -(interaction.elementIndices[1] - sample.get(0).get().size() + 1)
        //							: interaction.elementIndices[1] + 1;
        //					return array;
        //				}).filter(interaction -> {
        //
        //					for (BooleanSolution config : sample) {
        //						if (((interaction[0] < 0 && config.get(-interaction[0] - 1) == interaction[0])
        //								|| (interaction[0] > 0 && config.get(interaction[0] - 1) == interaction[0]))
        //								&& ((interaction[1] < 0 && config.get(-interaction[1] - 1) == interaction[1])
        //										|| (interaction[1] > 0 && config.get(interaction[1] - 1) == interaction[1]))) {
        //							return true;
        //						}
        //					}
        //					return false;
        //				}).count();
        //		System.out.println("Contained interactions in sample:" + countAppearingInteractions);
        //
        //		// count how many of the interactions in the sample also appear in the feature
        //		// model
        //		long countValidToFmSampleInteractions = LexicographicIterator.<Void>stream(t, sample.get(0).get().size() *
        // 2)
        //				.map(interaction -> {
        //					int[] array = new int[2];
        //					array[0] = interaction.elementIndices[0] > sample.get(0).get().size() - 1
        //							? -(interaction.elementIndices[0] - sample.get(0).get().size() + 1)
        //							: interaction.elementIndices[0] + 1;
        //					array[1] = interaction.elementIndices[1] > sample.get(0).get().size() - 1
        //							? -(interaction.elementIndices[1] - sample.get(0).get().size() + 1)
        //							: interaction.elementIndices[1] + 1;
        //					return array;
        //				}).filter(interaction -> {
        //
        //					SAT4JSolver solver = new SAT4JSolutionSolver(clauses.compute());
        //					VariableMap fMfeatures = Computations
        //							.await((IComputation<VariableMap>) cnf.map(Computations::getValue));
        //
        //					Result<String> fmFeatureName1 = features.get(Math.abs(interaction[0]));
        //					Result<String> fmFeatureName2 = features.get(Math.abs(interaction[1]));
        //
        //					if (!fMfeatures.get(fmFeatureName1.get()).isEmpty()
        //							&& !fMfeatures.get(fmFeatureName2.get()).isEmpty()) {
        //						for (BooleanSolution config : sample) {
        //							if (((interaction[0] < 0 && config.get(-interaction[0] - 1) == interaction[0])
        //									|| (interaction[0] > 0 && config.get(interaction[0] - 1) == interaction[0]))
        //									&& ((interaction[1] < 0 && config.get(-interaction[1] - 1) == interaction[1])
        //											|| (interaction[1] > 0
        //													&& config.get(interaction[1] - 1) == interaction[1]))) {
        //								int[] featureInteractions = new int[2];
        //								featureInteractions[0] = interaction[0] > 0 ? fMfeatures.get(fmFeatureName1.get()).get()
        //										: -fMfeatures.get(fmFeatureName1.get()).get();
        //								featureInteractions[1] = interaction[1] > 0 ? fMfeatures.get(fmFeatureName2.get()).get()
        //										: -fMfeatures.get(fmFeatureName2.get()).get();
        //
        //								if (solver.hasSolution(featureInteractions).get()) {
        //									return true;
        //								}
        //							}
        //						}
        //					}
        //
        //					return false;
        //				}).count();
        //
        //		System.out.println("Contained interactions in sample that are valid on fm:" +
        // countValidToFmSampleInteractions);

    }

    public static void calculateFMInteractions(int t, IComputation<BooleanClauseList> clauses) {

        // t = 1
        // count how many valid feature interactions exist overall in the feature model
        long countValidFMInteractions = LexicographicIterator.<Void>stream(
                        t, clauses.compute().getVariableCount() * 2)
                .map(interaction -> {
                    int[] array = new int[2];
                    array[0] = interaction.elementIndices[0] > clauses.compute().getVariableCount() - 1
                            ? -(interaction.elementIndices[0]
                                    - clauses.compute().getVariableCount()
                                    + 1)
                            : interaction.elementIndices[0] + 1;
                    return array;
                })
                .filter(interaction -> {
                    SAT4JSolver solver = new SAT4JSolutionSolver(clauses.compute());
                    if (solver.hasSolution(interaction).get()) {
                        return true;
                    }

                    return false;
                })
                .count();

        // t = 2
        // count how many valid feature interactions exist overall in the feature model
        //		long countValidFMInteractions = LexicographicIterator.<Void>stream(t, clauses.compute().getVariableCount() *
        // 2)
        //				.map(interaction -> {
        //					int[] array = new int[2];
        //					array[0] = interaction.elementIndices[0] > clauses.compute().getVariableCount() - 1
        //							? -(interaction.elementIndices[0] - clauses.compute().getVariableCount() + 1)
        //							: interaction.elementIndices[0] + 1;
        //					array[1] = interaction.elementIndices[1] > clauses.compute().getVariableCount() - 1
        //							? -(interaction.elementIndices[1] - clauses.compute().getVariableCount() + 1)
        //							: interaction.elementIndices[1] + 1;
        //					return array;
        //				}).filter(interaction -> {
        //
        //					SAT4JSolver solver = new SAT4JSolutionSolver(clauses.compute());
        //					if (solver.hasSolution(interaction).get()) {
        //						return true;
        //					}
        //
        //					return false;
        //				}).count();

        System.out.println("Valid interactions of fm:" + countValidFMInteractions);
    }
}
