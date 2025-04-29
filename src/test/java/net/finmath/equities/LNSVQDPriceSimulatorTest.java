package net.finmath.equities;

import net.finmath.equities.Simulation.Options.CliquetSimulationPricer;
import net.finmath.equities.Simulation.Options.EuropeanSimulationPricer;
import net.finmath.equities.Simulation.PathSimulator;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.modelling.descriptor.AssetModelDescriptor;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class LNSVQDPriceSimulatorTest extends TestsSetupForLNSVQD {
	Random random = new Random();

	/*@Test
	public void getImpliedVolsLNSVQD() throws Exception {
		// Set the right case
		loadS24(); //loadBitcoin(); //loadBitcoin(); // loadS24();

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size() / maturityGrid.length;

		// Get analytical prices
		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs);

		List<Integer> seeds = random.ints(10).boxed().collect(Collectors.toList());

		double[][][] pricesMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];
		double[][][] pricesQMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];

		double maxMaturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();

		int startingIndex = 0;
		double[] startingValueLNSVQD = new double[]{spot0, selectedParamsLNSVQD[0]};

		for(int seed : seeds) {
			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			pathSimulatorMC.precalculatePaths(seed, true, startingIndex + 1, startingValueLNSVQD, Boolean.TRUE);

			// Define pricers
			EuropeanSimulationPricer simulPricerMC = new EuropeanSimulationPricer<>(pathSimulatorMC);

			for(int m = 0; m < maturityGrid.length; m++) {
				double maturity = maturityGrid[m];
				double[] maturityGridQMC = new double[]{maturity};

				// QMC
				double[] timeGridQMC = LNSVQDUtils.addTimePointsToArray(maturityGridQMC,
								(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
						.stream().distinct().mapToDouble(Double::doubleValue).toArray();
				LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGridQMC, maturityGridQMC, lnsvqdModelAnalyticalPricer, false);
				pathSimulatorQMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD, Boolean.TRUE);

				// Define pricer
				EuropeanSimulationPricer simulPricerQMC = new EuropeanSimulationPricer<>(pathSimulatorQMC);

				for(int s = 0; s < numStrikesPerMaturity; s++) {
					double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();

					double simulatedOptionPrice;
					try {
						simulatedOptionPrice = simulPricerMC.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPrice = 1000000;
					}
					pricesMC[seeds.indexOf(seed)][m][s] = simulatedOptionPrice;

					// QMC
					double simulatedOptionPriceQMC;
					try {
						simulatedOptionPriceQMC = simulPricerQMC.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPriceQMC = 1000000;
					}
					pricesQMC[seeds.indexOf(seed)][m][s] = simulatedOptionPriceQMC;
				}
			}
			System.out.println(pricesMC[seeds.indexOf(seed)][4][3]);
			System.out.println("Finished seed " + seed);
		}

		for(int m = 0; m < maturityGrid.length; m++) {
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				double maturity = strikeMatPairs.get(m * numStrikesPerMaturity + s).getKey();
				double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();

				double[] pricesMCForPair = new double[seeds.size()];
				double[] pricesQMCForPair = new double[seeds.size()];
				for(int j = 0; j < seeds.size(); j++) {
					pricesMCForPair[j] = pricesMC[j][m][s];
					pricesQMCForPair[j] = pricesQMC[j][m][s];
				}

				double averagePrice = Math.max(Arrays.stream(pricesMCForPair).average().getAsDouble(), 1E-15);
				double varMC = Arrays.stream(pricesMCForPair).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
				double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalMC = LNSVQDUtils.getConfidenceInterval(pricesMCForPair, 0.05);
				double lbMc = Math.max(confidenceIntervalMC[0], 1E-10);
				double ubMc = Math.max(confidenceIntervalMC[1], 1E-10);

				double impliedVolMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePrice);
				double impliedVolLowerMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, lbMc);
				double impliedVolUpperMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, ubMc);

				double averagePriceQMC = Math.max(Arrays.stream(pricesQMCForPair).average().getAsDouble(), 1E-10);
				double varQMC = Arrays.stream(pricesQMCForPair).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
				double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalQMC = LNSVQDUtils.getConfidenceInterval(pricesQMCForPair, 0.05);
				double lbQmc = Math.max(confidenceIntervalQMC[0], 1E-10);
				double ubQmc = Math.max(confidenceIntervalQMC[1], 1E-10);

				double impliedVolQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePriceQMC);
				double impliedVolLowerQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, lbQmc);
				double impliedVolUpperQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, ubQmc);

				System.out.println(volAna[m * numStrikesPerMaturity + s] + "\t"
						+ impliedVolMC + "\t" + stdErrMC + "\t" + impliedVolLowerMC + "\t" + impliedVolUpperMC + "\t"
						+ impliedVolQMC + "\t" + stdErrQMC + "\t" + impliedVolLowerQMC + "\t" + impliedVolUpperQMC + "\t");

				*//*System.out.println(volAna[m * numStrikesPerMaturity + s] + "\t"
						+ averagePrice + "\t" + stdErrMC + "\t" + lbMc + "\t" + ubMc + "\t"
						+ averagePriceQMC + "\t" + stdErrQMC + "\t" + lbQmc + "\t" + ubQmc + "\t");*//*
			}
		}
	}*/

	@Test
	public void testCliquetOption() throws Exception {
		numberOfPaths = 100;

		// Set the right case
		loadS20(); // loadS20(); //loadS24(); //loadBitcoin(); // loadS24();

		// Set Cliquet params
		double maturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double floorL = Double.NEGATIVE_INFINITY;
		double capL = 1.1;
		double floorG = 0;
		double capG = Double.POSITIVE_INFINITY;

		List<Integer> seeds = random.ints(10).boxed().collect(Collectors.toList());

		double[] pricesLnsvqdMC = new double[seeds.size()];
		double[] pricesLnsvqdQMC = new double[seeds.size()];
		double[] pricesHestonMC = new double[seeds.size()];
		double[] pricesHestonQMC = new double[seeds.size()];

		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();

		int startingIndex = 0;
		double[] startingValueLNSVQD = new double[]{spot0, selectedParamsLNSVQD[0]};
		double[] startingValueHeston = new double[]{spot0, selectedParamsHeston[0]};

		for(int seed : seeds) {
			PathSimulator pathSimulatorLnsvqdMc = new PathSimulator(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid);
			PathSimulator pathSimulatorLnsvqdQmc = new PathSimulator(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid);
			PathSimulator pathSimulatorHestonMC = new PathSimulator(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid);
			PathSimulator pathSimulatorHestonQMC = new PathSimulator(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid);


			pathSimulatorLnsvqdMc.precalculatePaths(seed, true, startingIndex,
					startingValueLNSVQD, Boolean.TRUE, "LNSVQD", "MC", lnsvqdModelDescriptor, null);
			pathSimulatorLnsvqdQmc.precalculatePaths(seed, true, startingIndex,
					startingValueLNSVQD, Boolean.TRUE, "LNSVQD", "QMC", lnsvqdModelDescriptor, null);
			pathSimulatorHestonMC.precalculatePaths(seed, true, startingIndex,
					startingValueHeston, Boolean.TRUE, "Heston", "MC", null, hestonModelDescriptor);
			pathSimulatorHestonQMC.precalculatePaths(seed, true, startingIndex,
					startingValueHeston, Boolean.TRUE, "Heston", "QMC", null, hestonModelDescriptor);


			// Define pricers
			CliquetSimulationPricer simulPricerLnsvqdMC = new CliquetSimulationPricer(pathSimulatorLnsvqdMc);
			CliquetSimulationPricer simulPricerLnsvqdQMC = new CliquetSimulationPricer(pathSimulatorLnsvqdQmc);
			CliquetSimulationPricer simulPricerHestonMC = new CliquetSimulationPricer(pathSimulatorHestonMC);
			CliquetSimulationPricer simulPricerHestonQMC = new CliquetSimulationPricer(pathSimulatorHestonQMC);

			// MC
			double simulatedOptionPrice;
			try {
				simulatedOptionPrice = simulPricerLnsvqdMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPrice = 1000000;
			}
			pricesLnsvqdMC[seeds.indexOf(seed)] = simulatedOptionPrice;

			// QMC
			double simulatedOptionPriceQMC;
			try {
				simulatedOptionPriceQMC = simulPricerLnsvqdQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceQMC = 1000000;
			}
			pricesLnsvqdQMC[seeds.indexOf(seed)] = simulatedOptionPriceQMC;

			// MC Heston
			double simulatedOptionPriceHestonMC;
			try {
				simulatedOptionPriceHestonMC = simulPricerHestonMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonMC = 1000000;
			}
			pricesHestonMC[seeds.indexOf(seed)] = simulatedOptionPriceHestonMC;

			// QMC Heston
			double simulatedOptionPriceHestonQMC;
			try {
				simulatedOptionPriceHestonQMC = simulPricerHestonQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonQMC = 1000000;
			}
			pricesHestonQMC[seeds.indexOf(seed)] = simulatedOptionPriceHestonQMC;

			System.out.println(simulatedOptionPriceQMC);
			System.out.println("Finished seed " + seed);
		}

		double averagePrice = Math.max(Arrays.stream(pricesLnsvqdMC).average().getAsDouble(), 1E-10);
		double varMC = Arrays.stream(pricesLnsvqdMC).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
		double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());
		double[] confidenceIntervalMC = LNSVQDUtils.getConfidenceInterval(pricesLnsvqdMC, 0.05);
		double lbMC = Math.max(confidenceIntervalMC[0], 1E-10);
		double ubMC = Math.max(confidenceIntervalMC[1], 1E-10);

		double averagePriceQMC = Math.max(Arrays.stream(pricesLnsvqdQMC).average().getAsDouble(), 1E-10);
		double varQMC = Arrays.stream(pricesLnsvqdQMC).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
		double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());
		double[] confidenceIntervalQMC = LNSVQDUtils.getConfidenceInterval(pricesLnsvqdQMC, 0.05);
		double lbQMC = Math.max(confidenceIntervalQMC[0], 1E-10);
		double ubQMC = Math.max(confidenceIntervalQMC[1], 1E-10);

		double averagePriceHestonMC = Math.max(Arrays.stream(pricesHestonMC).average().getAsDouble(), 1E-10);
		double varHestonMC = Arrays.stream(pricesHestonMC).map(x -> Math.pow(x - averagePriceHestonMC, 2)).sum() / (seeds.size() - 1);
		double stdErrHestonMC = Math.sqrt(varHestonMC) / Math.sqrt(seeds.size());
		double[] confidenceIntervalHestonMC = LNSVQDUtils.getConfidenceInterval(pricesHestonMC, 0.05);
		double lbHestonMC = Math.max(confidenceIntervalHestonMC[0], 1E-10);
		double ubHestonMC = Math.max(confidenceIntervalHestonMC[1], 1E-10);

		double averagePriceHestonQMC = Math.max(Arrays.stream(pricesHestonQMC).average().getAsDouble(), 1E-10);
		double varHestonQMC = Arrays.stream(pricesHestonQMC).map(x -> Math.pow(x - averagePriceHestonQMC, 2)).sum() / (seeds.size() - 1);
		double stdErrHestonQMC = Math.sqrt(varHestonQMC) / Math.sqrt(seeds.size());
		double[] confidenceIntervalHestonQMC = LNSVQDUtils.getConfidenceInterval(pricesHestonQMC, 0.05);
		double lbHestonQMC = Math.max(confidenceIntervalHestonQMC[0], 1E-10);
		double ubHestonQMC = Math.max(confidenceIntervalHestonQMC[1], 1E-10);

		System.out.println(averagePrice + "\t" + stdErrMC + "\t" + lbMC + "\t" + ubMC + "\t"
				+ averagePriceQMC + "\t" + stdErrQMC + "\t" + lbQMC + "\t" + ubQMC + "\t"
				+ averagePriceHestonMC + "\t" + stdErrHestonMC + "\t" + lbHestonMC + "\t" + ubHestonMC + "\t"
				+ averagePriceHestonQMC + "\t" + stdErrHestonQMC + "\t" + lbHestonQMC + "\t" + ubHestonQMC);
	}

}