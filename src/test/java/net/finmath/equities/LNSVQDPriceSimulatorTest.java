package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.*;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import java.util.ArrayList;
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
		ArrayList<Pair<Double, Double>> strikeMatPairs = setDAXHestonSetupSIM(); //setBTCSetupSIM(); // setDAXHestonSetupSIM();

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size() / maturityGrid.length;

		// Get analytical prices
		StopWatch sw = StopWatch.createStarted();
		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs, null);
		sw.stop();
		System.out.println("time: " + sw.getTime());

		List<Integer> seeds = random.ints(5).boxed().collect(Collectors.toList());

		double[][][] pricesMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];
		double[][][] pricesQMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];

		double maxMaturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		for(int seed : seeds) {
			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, maturityGrid, false);
			// sw.reset();
			// sw.start();
			pathSimulatorMC.precalculatePaths(seed, true);
			// System.out.println("time MC: " + sw.getTime());

			// Define pricers
			LNSVQDEuropeanSimulationPricer simulPricerMC = new LNSVQDEuropeanSimulationPricer<>(pathSimulatorMC);

			for(int m = 0; m < maturityGrid.length; m++) {
				double maturity = maturityGrid[m];
				double[] maturityGridQMC = new double[]{maturity};

				// QMC
				double[] timeGridQMC = LNSVQDUtils.addTimePointsToArray(maturityGridQMC,
								(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
						.stream().distinct().mapToDouble(Double::doubleValue).toArray();
				LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGridQMC, maturityGridQMC, false);
				// sw.reset();
				// sw.start();
				pathSimulatorQMC.precalculatePaths(seed, true);
				// sw.stop();
				// System.out.println("time QMC: " + sw.getTime());

				// Define pricer
				LNSVQDEuropeanSimulationPricer simulPricerQMC = new LNSVQDEuropeanSimulationPricer<>(pathSimulatorQMC);

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

				double averagePrice = Math.max(Arrays.stream(pricesMCForPair).average().getAsDouble(), 1E-10);
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
			}
		}
	}

	@Test
	public void testCliquetOption() throws Exception {
		// Set the right case
		ArrayList<Pair<Double, Double>> strikeMatPairs = setDAXHestonSetupSIM(); //setBTCSetupSIM(); // setDAXHestonSetupSIM();

		// Set Cliquet params
		double maturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double floorL = 0;
		double capL = 10;
		double floorG = 0;
		double capG = 10;

		List<Integer> seeds = random.ints(2).boxed().collect(Collectors.toList());

		double[] pricesMC = new double[seeds.size()];
		double[] pricesQMC = new double[seeds.size()];

		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		for(int seed : seeds) {
			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, maturityGrid, false);
			// sw.reset();
			// sw.start();
			pathSimulatorMC.precalculatePaths(seed, true);
			// System.out.println("time MC: " + sw.getTime());

			LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, maturityGrid, false);
			// sw.reset();
			// sw.start();
			pathSimulatorQMC.precalculatePaths(seed, true);
			// sw.stop();
			// System.out.println("time QMC: " + sw.getTime());

			// Define pricers
			LNSVQDCliquetPricer simulPricerMC = new LNSVQDCliquetPricer<>(pathSimulatorMC);
			LNSVQDCliquetPricer simulPricerQMC = new LNSVQDCliquetPricer<>(pathSimulatorQMC);

			double simulatedOptionPrice;
			try {
				simulatedOptionPrice = simulPricerMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPrice = 1000000;
			}
			pricesMC[seeds.indexOf(seed)] = simulatedOptionPrice;

			// QMC
			double simulatedOptionPriceQMC;
			try {
				simulatedOptionPriceQMC = simulPricerQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceQMC = 1000000;
			}
			pricesQMC[seeds.indexOf(seed)] = simulatedOptionPriceQMC;

			System.out.println("Finished seed " + seed);
		}

		double averagePrice = Math.max(Arrays.stream(pricesMC).average().getAsDouble(), 1E-10);
		double varMC = Arrays.stream(pricesMC).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
		double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());
		double[] confidenceIntervalMC = LNSVQDUtils.getConfidenceInterval(pricesMC, 0.05);
		double lbMc = Math.max(confidenceIntervalMC[0], 1E-10);
		double ubMc = Math.max(confidenceIntervalMC[1], 1E-10);

		double averagePriceQMC = Math.max(Arrays.stream(pricesQMC).average().getAsDouble(), 1E-10);
		double varQMC = Arrays.stream(pricesQMC).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
		double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());
		double[] confidenceIntervalQMC = LNSVQDUtils.getConfidenceInterval(pricesQMC, 0.05);
		double lbQmc = Math.max(confidenceIntervalQMC[0], 1E-10);
		double ubQmc = Math.max(confidenceIntervalQMC[1], 1E-10);

		System.out.println(averagePrice + "\t" + stdErrMC + "\t" + lbMc + "\t" + ubMc + "\t"
				+ averagePriceQMC + "\t" + stdErrQMC + "\t" + lbQmc + "\t" + ubQmc);
	}

}