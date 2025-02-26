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
	public void testCase() throws Exception {
		// Set the right case
		ArrayList<Pair<Double, Double>> strikeMatPairs = setDAXHestonSetupSIM(); //setBTCSetupSIM(); // setDAXHestonSetupSIM();

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size() / maturityGrid.length;

		// Get analytical prices
		StopWatch sw = StopWatch.createStarted();
		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs, null);
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
			LNSVQDEuropeanPriceSimulator lnsvqdPriceSimulator = new LNSVQDEuropeanPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, maturityGrid, false);
			// sw.reset();
			// sw.start();
			lnsvqdPriceSimulator.precalculatePaths(seed, true);
			// System.out.println("time MC: " + sw.getTime());

			LNSVQDEuropeanPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDEuropeanPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, maturityGrid, false);
			// sw.reset();
			// sw.start();
			lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
			// sw.stop();
			// System.out.println("time MC: " + sw.getTime());

			for(int m = 0; m < maturityGrid.length; m++) {
				for(int s = 0; s < numStrikesPerMaturity; s++) {
					double maturity = strikeMatPairs.get(m * numStrikesPerMaturity + s).getKey();
					double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();

					double simulatedOptionPrice;
					try {
						simulatedOptionPrice = lnsvqdPriceSimulator.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPrice = 1000000;
					}
					pricesMC[seeds.indexOf(seed)][m][s] = simulatedOptionPrice;

					// QMC
					double simulatedOptionPriceQMC;
					try {
						simulatedOptionPriceQMC = lnsvqdPriceSimulatorQMC.getEuropeanPriceAuto(strike, maturity);
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
	public void testDependencyOnDiscretization() throws Exception {
		// Set the right case
		ArrayList<Pair<Double, Double>> strikeMatPairs = setDAXHestonSetupSIM(); //setBTCSetupSIM(); // setDAXHestonSetupSIM();

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size() / maturityGrid.length;

		List<Integer> seeds = IntStream.range(1, 6).boxed().collect(Collectors.toList());

		double[][][] pricesMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];
		double[][][] pricesQMC = new double[seeds.size()][maturityGrid.length][numStrikesPerMaturity];

		for(int seed : seeds) {
			double maxMaturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
			double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
							(int) (Math.round(maxMaturity * 365.) * seed), 0, maxMaturity, true)
					.stream().distinct().mapToDouble(Double::doubleValue).toArray();

			LNSVQDEuropeanPriceSimulator lnsvqdPriceSimulator = new LNSVQDEuropeanPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, maturityGrid, false);
			lnsvqdPriceSimulator.precalculatePaths(5603, true);

			LNSVQDEuropeanPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDEuropeanPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, maturityGrid, false);
			// lnsvqdPriceSimulatorQMC.precalculatePaths(5603);

			for(int m = 0; m < maturityGrid.length; m++) {
				for(int s = 0; s < numStrikesPerMaturity; s++) {
					double maturity = strikeMatPairs.get(m * numStrikesPerMaturity + s).getKey();
					double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();

					double simulatedOptionPrice;
					try {
						simulatedOptionPrice = lnsvqdPriceSimulator.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPrice = 1000000;
					}
					pricesMC[seeds.indexOf(seed)][m][s] = simulatedOptionPrice;

					// QMC
					double simulatedOptionPriceQMC;
					try {
						simulatedOptionPriceQMC = 0; // lnsvqdPriceSimulatorQMC.getEuropeanPriceAuto(strike, maturity);
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
					pricesMCForPair[j] = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, pricesMC[j][m][s]);
					pricesQMCForPair[j] = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, pricesQMC[j][m][s]);
				}
				LNSVQDUtils.printArray(pricesMCForPair);
				LNSVQDUtils.printArray(pricesQMCForPair);
			}
		}
	}

	/*@Test
	public void printTransformedVolPathMoments() throws CalculationException {
		int numberOfPaths = 100000;
		// Get option values
		double maturity = 1.2;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
				maturity, (int) Math.round(maturity * 365.));

		int seed = 1;

		LNSVQDEuropeanPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDEuropeanPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
		lnsvqdCallPriceSimulator.precalculatePaths(seed);
		double[][][] transformedPaths = lnsvqdCallPriceSimulator.getTransformedPath();

		for(int j = 0; j < timeGrid.length; j++) {
			double moment1 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 1)).average().getAsDouble();
			double moment2 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 2)).average().getAsDouble();
			double moment3 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 3)).average().getAsDouble();
			double moment4 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 4)).average().getAsDouble();
			double moment5 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 5)).average().getAsDouble();
			double moment6 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 6)).average().getAsDouble();
			System.out.println(timeGrid[j] +
					"\t" + moment1 +
					"\t" + moment2 +
					"\t" + moment3 +
					"\t" + moment4 +
					"\t" + moment5 +
					"\t" + moment6);
		}
	}*/
}