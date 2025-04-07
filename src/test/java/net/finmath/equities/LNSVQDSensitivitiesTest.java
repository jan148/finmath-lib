package net.finmath.equities;

import net.finmath.equities.Simulation.HestonPathSimulator.HestonPathSimulatorMC;
import net.finmath.equities.Simulation.HestonPathSimulator.HestonPathSimulatorQMC;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorMC;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorQMC;
import net.finmath.equities.Simulation.Options.CliquetSimulationPricer;
import net.finmath.equities.Simulation.Options.EuropeanSimulationPricer;
import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.LNSVQDUtils;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class LNSVQDSensitivitiesTest extends TestsSetupForLNSVQD {
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
		StopWatch sw = StopWatch.createStarted();
		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs);
		sw.stop();
		System.out.println("time: " + sw.getTime());

		List<Integer> seeds = random.ints(10).boxed().collect(Collectors.toList());

		double[][][] pricesMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];
		double[][][] pricesQMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];

		double maxMaturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		for(int seed : seeds) {
			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			// sw.reset();
			// sw.start();
			pathSimulatorMC.precalculatePaths(seed, true);
			// System.out.println("time MC: " + sw.getTime());

			// Define pricers
			EuropeanSimulationPricer simulPricerMC = new EuropeanSimulationPricer<>(pathSimulatorMC);

			for(int m = 0; m < maturityGrid.length; m++) {
				double maturity = maturityGrid[m];
				double[] maturityGridQMC = new double[]{maturity};

				// QMC
				double[] timeGridQMC = LNSVQDUtils.addTimePointsToArray(maturityGridQMC,
								(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
						.stream().distinct().mapToDouble(Double::doubleValue).toArray();
				LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
				// sw.reset();
				// sw.start();
				// pathSimulatorQMC.precalculatePaths(seed, true);
				// sw.stop();
				// System.out.println("time QMC: " + sw.getTime());

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
	public void testCliquetDelta() throws Exception {
		// Set the right case
		setDAXHestonFebruarySetupSIM();

		// Set Cliquet params
		double maturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double floorL = Double.NEGATIVE_INFINITY;
		double capL = 1.1;
		double floorG = 0;
		double capG = Double.POSITIVE_INFINITY;

		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();

		int minPercOfSpot = 95;
		int maxPercOfSpot = 105;
		double shiftSize = 1E-8;
		double[] S = IntStream.range(minPercOfSpot, maxPercOfSpot + 1).mapToDouble(n -> spot0 * (n / 100.)).toArray();
		double[] sShifted = Arrays.stream(S).map(s -> s + shiftSize).toArray();

		double[] pricesLnsvqdMC = new double[S.length];
		double[] pricesLnsvqdQMC = new double[S.length];
		double[] pricesHestonMC = new double[S.length];
		double[] pricesHestonQMC = new double[S.length];

		double[] pricesLnsvqdMCShifted = new double[S.length];
		double[] pricesLnsvqdQMCShifted = new double[S.length];
		double[] pricesHestonMCShifted = new double[S.length];
		double[] pricesHestonQMCShifted= new double[S.length];

		for(int k = 0; k < S.length; k++) {
			int seed = 2956;
			double spot = S[k];
			lnsvqdModelAnalyticalPricer.setSpot(spot);
			equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, spot, forwardCurve, affineDividendStream, dayCountConvention);

			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			// sw.reset();
			// sw.start();
			pathSimulatorMC.precalculatePaths(seed, true);
			// System.out.println("time MC: " + sw.getTime());

			// QMC
			LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			// sw.reset();
			// sw.start();
			pathSimulatorQMC.precalculatePaths(seed, true);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Heston MC
			double gamma1 = 1;
			double gamma2 = 0;
			HestonPathSimulatorMC pathSimulatorHestonMC = new HestonPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
			// sw.reset();
			// sw.start();
			pathSimulatorHestonMC.precalculatePaths(seed, true);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Heston QMC
			HestonPathSimulatorQMC pathSimulatorHestonQMC = new HestonPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);

			// sw.reset();
			// sw.start();
			pathSimulatorHestonQMC.precalculatePaths(seed, true);
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
			pricesLnsvqdMC[k] = simulatedOptionPrice;

			// QMC
			double simulatedOptionPriceQMC;
			try {
				simulatedOptionPriceQMC = simulPricerLnsvqdQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceQMC = 1000000;
			}
			pricesLnsvqdQMC[k] = simulatedOptionPriceQMC;

			// MC Heston
			double simulatedOptionPriceHestonMC;
			try {
				simulatedOptionPriceHestonMC = simulPricerHestonMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonMC = 1000000;
			}
			pricesHestonMC[k] = simulatedOptionPriceHestonMC;

			// QMC Heston
			double simulatedOptionPriceHestonQMC;
			try {
				simulatedOptionPriceHestonQMC = simulPricerHestonQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonQMC = 1000000;
			}
			pricesHestonQMC[k] = simulatedOptionPriceHestonQMC;

			System.out.println("Finished spot " + spot);
		}

		for(int k = 0; k < sShifted.length; k++) {
			int seed = 2956;
			double spot = sShifted[k];
			lnsvqdModelAnalyticalPricer.setSpot(spot);
			equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, spot, forwardCurve, affineDividendStream, dayCountConvention);

			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);;
			// sw.reset();
			// sw.start();
			pathSimulatorMC.precalculatePaths(seed, true);
			// System.out.println("time MC: " + sw.getTime());

			// QMC
			LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			// sw.reset();
			// sw.start();
			pathSimulatorQMC.precalculatePaths(seed, true);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Heston MC
			double gamma1 = 1;
			double gamma2 = 0;
			HestonPathSimulatorMC pathSimulatorHestonMC = new HestonPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
			// sw.reset();
			// sw.start();
			pathSimulatorHestonMC.precalculatePaths(seed, true);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Heston QMC
			HestonPathSimulatorQMC pathSimulatorHestonQMC = new HestonPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);

			// sw.reset();
			// sw.start();
			pathSimulatorHestonQMC.precalculatePaths(seed, true);
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
			pricesLnsvqdMCShifted[k] = simulatedOptionPrice;

			// QMC
			double simulatedOptionPriceQMC;
			try {
				simulatedOptionPriceQMC = simulPricerLnsvqdQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceQMC = 1000000;
			}
			pricesLnsvqdQMCShifted[k] = simulatedOptionPriceQMC;

			// MC Heston
			double simulatedOptionPriceHestonMC;
			try {
				simulatedOptionPriceHestonMC = simulPricerHestonMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonMC = 1000000;
			}
			pricesHestonMCShifted[k] = simulatedOptionPriceHestonMC;

			// QMC Heston
			double simulatedOptionPriceHestonQMC;
			try {
				simulatedOptionPriceHestonQMC = simulPricerHestonQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonQMC = 1000000;
			}
			pricesHestonQMCShifted[k] = simulatedOptionPriceHestonQMC;

			System.out.println("Finished spot " + spot);
		}

		double[] firstOrderMC = getSensitivities(S, pricesLnsvqdMC, sShifted, pricesLnsvqdMCShifted);
		double[] firstOrderQMC = getSensitivities(S, pricesLnsvqdQMC, sShifted, pricesLnsvqdQMCShifted);
		double[] firstOrderHestonMC = getSensitivities(S, pricesHestonMC, sShifted, pricesHestonMCShifted);
		double[] firstOrderHestonQMC = getSensitivities(S, pricesHestonQMC, sShifted, pricesHestonQMCShifted);

		LNSVQDUtils.printArray(S);
		LNSVQDUtils.printArray(firstOrderMC);
		LNSVQDUtils.printArray(firstOrderQMC);
		LNSVQDUtils.printArray(firstOrderHestonMC);
		LNSVQDUtils.printArray(firstOrderHestonQMC);
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
			pathSimulatorMC.precalculatePaths(seed, true);
			// System.out.println("time MC: " + sw.getTime());

			// Heston MC
			double gamma1 = 1;
			double gamma2 = 0;
			HestonPathSimulatorMC pathSimulatorHestonMC = new HestonPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
			// sw.reset();
			// sw.start();
			pathSimulatorHestonMC.precalculatePaths(seed, true);
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

	private double[] getSensitivities(double[] X, double[] Y, double[] xShifted, double[] yShifted) {
		double[] firstOrderDerivatives = new double[X.length];
		for(int k = 0; k < X.length; k++) {
			firstOrderDerivatives[k] = (yShifted[k] - Y[k]) / (xShifted[k] - X[k]);
		}
		return firstOrderDerivatives;
	}
}