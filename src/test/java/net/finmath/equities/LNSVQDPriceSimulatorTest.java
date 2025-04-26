package net.finmath.equities;

import net.finmath.equities.Simulation.HestonPathSimulator.HestonPathSimulatorMC;
import net.finmath.equities.Simulation.HestonPathSimulator.HestonPathSimulatorQMC;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorMC;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorQMC;
import net.finmath.equities.Simulation.Options.CliquetSimulationPricer;
import net.finmath.equities.Simulation.Options.EuropeanSimulationPricer;
import net.finmath.equities.models.LNSVQDUtils;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class LNSVQDPriceSimulatorTest extends TestsSetupForLNSVQD {
	/**
	 * Stat utils
	 */
	StandardDeviation standardDeviation = new StandardDeviation();
	Random random = new Random();

	/**
	 * Test cases: DAX: setDAXHestonSetupSIM(); BTC: setBTCSetupSIM()
	 */
	@Test
	public void testEuropeanOption() throws Exception {
		// Set the right case
		setDAXHestonSetupSIM(); //setBTCSetupSIM(); //setBTCSetupSIM(); // setDAXHestonSetupSIM();

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

		int startingIndex = 1;
		double[] startingValueLNSVQD = new double[]{spot0, selectedParamsLNSVQD[0]};
		double[] startingValueHeston = new double[]{spot0, selectedParamsHeston[0]};

		for(int seed : seeds) {
			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			pathSimulatorMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD, Boolean.TRUE);

			// Define pricers
			EuropeanSimulationPricer simulPricerMC = new EuropeanSimulationPricer<>(pathSimulatorMC);

			for(int m = 0; m < maturityGrid.length; m++) {
				double maturity = maturityGrid[m];
				double[] maturityGridQMC = new double[]{maturity};

				// QMC
				double[] timeGridQMC = LNSVQDUtils.addTimePointsToArray(maturityGridQMC,
								(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
						.stream().distinct().mapToDouble(Double::doubleValue).toArray();
				LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGridQMC, maturityGrid, lnsvqdModelAnalyticalPricer, false);
				// pathSimulatorQMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD);

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
						simulatedOptionPriceQMC = 0; // simulPricerQMC.getEuropeanPriceAuto(strike, maturity);
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

				/*System.out.println(volAna[m * numStrikesPerMaturity + s] + "\t"
						+ averagePrice + "\t" + stdErrMC + "\t" + lbMc + "\t" + ubMc + "\t"
						+ averagePriceQMC + "\t" + stdErrQMC + "\t" + lbQmc + "\t" + ubQmc + "\t");*/
			}
		}
	}

	@Test
	public void testCliquetOption() throws Exception {
		numberOfPaths = 100;

		// Set the right case
		setDAXHestonFebruarySetupSIM(); // setDAXHestonMarchSetupSIM(); //setDAXHestonSetupSIM(); //setBTCSetupSIM(); // setDAXHestonSetupSIM();

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

		int startingIndex = 1;
		double[] startingValueLNSVQD = new double[]{spot0, selectedParamsLNSVQD[0]};
		double[] startingValueHeston = new double[]{spot0, selectedParamsHeston[0]};

		for(int seed : seeds) {
			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);;
			// sw.reset();
			// sw.start();
			pathSimulatorMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD, Boolean.TRUE);
			// System.out.println("time MC: " + sw.getTime());

			// QMC
			LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			// sw.reset();
			// sw.start();
			pathSimulatorQMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD, Boolean.TRUE);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Heston MC
			double gamma1 = 1;
			double gamma2 = 0;
			HestonPathSimulatorMC pathSimulatorHestonMC = new HestonPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
			// sw.reset();
			// sw.start();
			pathSimulatorHestonMC.precalculatePaths(seed, true, startingIndex, startingValueHeston, Boolean.TRUE);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Heston QMC
			HestonPathSimulatorQMC pathSimulatorHestonQMC = new HestonPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);

			// sw.reset();
			// sw.start();
			pathSimulatorHestonQMC.precalculatePaths(seed, true, startingIndex, startingValueHeston, Boolean.TRUE);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Define pricers
			CliquetSimulationPricer simulPricerLnsvqdMC = new CliquetSimulationPricer<>(pathSimulatorMC);
			CliquetSimulationPricer simulPricerLnsvqdQMC = new CliquetSimulationPricer<>(pathSimulatorQMC);
			CliquetSimulationPricer simulPricerHestonMC = new CliquetSimulationPricer<>(pathSimulatorHestonMC);
			CliquetSimulationPricer simulPricerHestonQMC = new CliquetSimulationPricer<>(pathSimulatorHestonQMC);

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

	@Test
	public void testHestonEuropean() throws Exception {
		// Set the right case
		setDAXHestonSetupSIM(); //setBTCSetupSIM(); // setDAXHestonSetupSIM();

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size() / maturityGrid.length;

		// Get analytical prices
		StopWatch sw = StopWatch.createStarted();
		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs);
		double[] pricesAna = lnsvqdModelAnalyticalPricer.getEuropeanOptionPricesAuto(strikeMatPairs);
		sw.stop();
		System.out.println("time: " + sw.getTime());

		List<Integer> seeds = random.ints(10).boxed().collect(Collectors.toList());

		double[][][] pricesMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];
		double[][][] pricesHestonMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];

		double maxMaturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		for(int seed : seeds) {
			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			// sw.reset();
			// sw.start();
			pathSimulatorMC.precalculatePaths(seed, true, 0, null, Boolean.TRUE);
			// System.out.println("time MC: " + sw.getTime());

			// Heston MC
			double gamma1 = 1;
			double gamma2 = 0;
			HestonPathSimulatorMC pathSimulatorHestonMC = new HestonPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
			// sw.reset();
			// sw.start();
			pathSimulatorHestonMC.precalculatePaths(seed, true, 0, null, Boolean.TRUE);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Define pricers
			EuropeanSimulationPricer simulPricerMC = new EuropeanSimulationPricer<>(pathSimulatorMC);
			EuropeanSimulationPricer simulPricerHestonMC = new EuropeanSimulationPricer<>(pathSimulatorHestonMC);

			for(int m = 0; m < maturityGrid.length; m++) {
				double maturity = maturityGrid[m];

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

					// MC Heston
					double simulatedOptionPriceHestonMC;
					try {
						simulatedOptionPriceHestonMC = simulPricerHestonMC.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPriceHestonMC = 1000000;
					}
					pricesHestonMC[seeds.indexOf(seed)][m][s] = simulatedOptionPriceHestonMC;
				}
			}
			System.out.println("Finished seed " + seed);
		}

		for(int m = 0; m < maturityGrid.length; m++) {
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				double maturity = strikeMatPairs.get(m * numStrikesPerMaturity + s).getKey();
				double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();

				double[] pricesMCForPair = new double[seeds.size()];
				double[] pricesHestonMCForPair = new double[seeds.size()];
				for(int j = 0; j < seeds.size(); j++) {
					pricesMCForPair[j] = pricesMC[j][m][s];
					pricesHestonMCForPair[j] = pricesHestonMC[j][m][s];
				}

				double averagePrice = Math.max(Arrays.stream(pricesMCForPair).average().getAsDouble(), 1E-10);
				double varMC = Arrays.stream(pricesMCForPair).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
				double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalMC = LNSVQDUtils.getConfidenceInterval(pricesMCForPair, 0.05);
				double lbMc = Math.max(confidenceIntervalMC[0], 1E-10);
				double ubMc = Math.max(confidenceIntervalMC[1], 1E-10);

				double impliedVolMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePrice);
				double impliedVolLowerMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, lbMc);
				double impliedVolUpperMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, ubMc);

				double averagePriceHestonMC = Math.max(Arrays.stream(pricesHestonMCForPair).average().getAsDouble(), 1E-10);
				double varHestonMC = Arrays.stream(pricesHestonMCForPair).map(x -> Math.pow(x - averagePriceHestonMC, 2)).sum() / (seeds.size() - 1);
				double stdErrHestonMC = Math.sqrt(varHestonMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalHestonMC = LNSVQDUtils.getConfidenceInterval(pricesHestonMCForPair, 0.05);
				double lbHestonMc = Math.max(confidenceIntervalHestonMC[0], 1E-10);
				double ubHestonMc = Math.max(confidenceIntervalHestonMC[1], 1E-10);

				double impliedVolHestonMc = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePriceHestonMC);
				double impliedVolLowerHestonMc = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, lbHestonMc);
				double impliedVolUpperHestonMc = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, ubHestonMc);

				/*System.out.println(volAna[m * numStrikesPerMaturity + s] + "\t"
						+ impliedVolMC + "\t" + stdErrMC + "\t" + impliedVolLowerMC + "\t" + impliedVolUpperMC + "\t"
						+ impliedVolHestonMc + "\t" + stdErrHestonMC + "\t" + impliedVolLowerHestonMc + "\t" + impliedVolUpperHestonMc + "\t");
*/
				System.out.println(pricesAna[m * numStrikesPerMaturity + s] + "\t"
						+ averagePrice + "\t" + stdErrMC + "\t" + lbMc + "\t" + ubMc + "\t"
						+ averagePriceHestonMC + "\t" + stdErrHestonMC + "\t" + lbHestonMc + "\t" + ubHestonMc + "\t");
			}
		}
	}

}