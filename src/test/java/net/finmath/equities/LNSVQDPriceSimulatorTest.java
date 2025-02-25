package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.*;
import net.finmath.exception.CalculationException;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class LNSVQDPriceSimulatorTest extends TestsSetupForLNSVQD {
	/**
	 * Stat utils
	 */
	StandardDeviation standardDeviation = new StandardDeviation();
	Random random = new Random();

	@Test
	public void testDax() throws Exception {
		ArrayList<Pair<Double, Double>> strikeMatPairs = setDAXHestonSetup();

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size() / maturityGrid.length;

		double[][] pricesMC = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] pricesQMC = new double[maturityGrid.length][numStrikesPerMaturity];

		// Get analytical prices
		StopWatch sw = StopWatch.createStarted();
		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs, null);
		sw.stop();
		System.out.println("time: " + sw.getTime()); // formatted string like "12.3 ms"

		double[][] stdErrorsMc = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] stdErrorsQMc = new double[maturityGrid.length][numStrikesPerMaturity];

		for(int m = 0; m < maturityGrid.length; m++) {
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				double maturity = strikeMatPairs.get(m * numStrikesPerMaturity + s).getKey();
				double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();
				double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
						maturity, (int) (Math.round(maturity * 365.) * 2));

				List<Integer> seeds = random.ints(30).boxed().collect(Collectors.toList());
				double[] prices = new double[seeds.size()];
				double[] pricesQ = new double[seeds.size()];

				for(int seed : seeds) {
					// Normal MC
					// sw.reset();
					// sw.start();
					LNSVQDEuropeanPriceSimulator lnsvqdPriceSimulator = new LNSVQDEuropeanPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdPriceSimulator.precalculatePaths(seed);
					// sw.stop();
					// System.out.println("time MC: " + sw.getTime());
					double simulatedOptionPrice;
					try {
						simulatedOptionPrice = lnsvqdPriceSimulator.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPrice = 1000000;
					}

					prices[seeds.indexOf(seed)] = simulatedOptionPrice;

					// QMC
					LNSVQDPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					sw.reset();
					sw.start();
					// lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
					sw.stop();
					System.out.println("time MC: " + sw.getTime());
					double simulatedOptionPriceQMC;
					try {
						// simulatedOptionPriceQMC = lnsvqdPriceSimulatorQMC.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						// simulatedOptionPriceQMC = 1000000;
					}
					pricesQ[seeds.indexOf(seed)] = 0; // simulatedOptionPriceQMC;
				}

				double averagePrice = Arrays.stream(prices).average().getAsDouble();
				double varMC = Arrays.stream(prices).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
				double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());
				double under = averagePrice + stdErrMC * 1.96; //
				System.out.println(lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, under)); //
				double[] confidenceIntervalMC = LNSVQDUtils.getConfidenceInterval(prices, 0.05);

				double impliedVolMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePrice);
				double impliedVolLowerMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalMC[0]);
				double impliedVolUpperMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalMC[1]);

				double averagePriceQMC = Arrays.stream(pricesQ).average().getAsDouble();
				double varQMC = Arrays.stream(pricesQ).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
				double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalQMC = LNSVQDUtils.getConfidenceInterval(pricesQ, 0.05);

				double impliedVolQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePriceQMC);
				double impliedVolLowerQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalQMC[0]);
				double impliedVolUpperQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalQMC[1]);

				pricesMC[m][s] = averagePrice;
				pricesQMC[m][s] = averagePriceQMC;

				stdErrorsMc[m][s] = stdErrMC;
				stdErrorsQMc[m][s] = stdErrQMC;

				System.out.println(volAna[m * numStrikesPerMaturity + s] + "\t"
						+ impliedVolMC + "\t" + stdErrMC + "\t" + impliedVolLowerMC + "\t" + impliedVolUpperMC + "\t"
						+ impliedVolQMC  + "\t" + stdErrQMC + "\t" + impliedVolLowerQMC + "\t" + impliedVolUpperQMC + "\t");

				/*System.out.println(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t"
						+ pricesMC[m][s] + "\t" + stdErrMC + "\t" + confidenceIntervalMC[0] + "\t" + confidenceIntervalMC[1] + "\t"
						+ pricesQMC[m][s] + "\t" + stdErrQMC + "\t" + confidenceIntervalQMC[0] + "\t" + confidenceIntervalQMC[1] + "\t");*/
			}
		}
	}

	// Nest method assumes the same number of strikes for all maturities!
	@Test
	public void testBTC() throws Exception {
		ArrayList<Pair<Double, Double>> strikeMatPairs = setBTCSetup();

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size() / maturityGrid.length;

		double[][] pricesMC = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] pricesQMC = new double[maturityGrid.length][numStrikesPerMaturity];

		// Get analytical prices
		// StopWatch sw = StopWatch.createStarted();
		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs, null);
		// sw.stop();
		// System.out.println("time: " + sw.getTime()); // formatted string like "12.3 ms"

		double[][] stdErrorsMc = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] stdErrorsQMc = new double[maturityGrid.length][numStrikesPerMaturity];

		for(int m = 0; m < maturityGrid.length; m++) {
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				double maturity = strikeMatPairs.get(m * numStrikesPerMaturity + s).getKey();
				double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();
				double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
						maturity, (int) (Math.round(maturity * 365.) * 1));

				List<Integer> seeds = random.ints(10).boxed().collect(Collectors.toList());
				double[] prices = new double[seeds.size()];
				double[] pricesQ = new double[seeds.size()];

				for(int seed : seeds) {
					// Normal MC
					// sw.reset();
					// sw.start();
					LNSVQDEuropeanPriceSimulator lnsvqdPriceSimulator = new LNSVQDEuropeanPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdPriceSimulator.precalculatePaths(seed);
					// sw.stop();
					// System.out.println("time MC: " + sw.getTime());
					double simulatedOptionPrice;
					try {
						simulatedOptionPrice = lnsvqdPriceSimulator.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPrice = 1000000;
					}

					prices[seeds.indexOf(seed)] = simulatedOptionPrice;

					// QMC
					LNSVQDPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					// sw.reset();
					// sw.start();
					lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
					// sw.stop();
					// System.out.println("time MC: " + sw.getTime());
					double simulatedOptionPriceQMC;
					try {
						simulatedOptionPriceQMC = lnsvqdPriceSimulatorQMC.getEuropeanPriceAuto(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPriceQMC = 1000000;
					}
					pricesQ[seeds.indexOf(seed)] = simulatedOptionPriceQMC;
				}

				double averagePrice = Arrays.stream(prices).average().getAsDouble();
				double varMC = Arrays.stream(prices).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
				double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalMC = LNSVQDUtils.getConfidenceInterval(prices, 0.05);

				double impliedVolMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePrice);
				double impliedVolLowerMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalMC[0]);
				double impliedVolUpperMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalMC[1]);

				double averagePriceQMC = Arrays.stream(pricesQ).average().getAsDouble();
				double varQMC = Arrays.stream(pricesQ).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
				double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalQMC = LNSVQDUtils.getConfidenceInterval(pricesQ, 0.05);

				double impliedVolQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, averagePriceQMC);
				double impliedVolLowerQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalQMC[0]);
				double impliedVolUpperQMC = lnsvqdModelAnalyticalPricer.getImpliedVolFromPrice(strike, maturity, confidenceIntervalQMC[1]);

				pricesMC[m][s] = averagePrice;
				pricesQMC[m][s] = averagePriceQMC;

				stdErrorsMc[m][s] = stdErrMC;
				stdErrorsQMc[m][s] = stdErrQMC;

				System.out.println(volAna[m * numStrikesPerMaturity + s] + "\t"
						+ impliedVolMC + "\t" + stdErrMC + "\t" + impliedVolLowerMC + "\t" + impliedVolUpperMC + "\t"
						+ impliedVolQMC  + "\t" + stdErrQMC + "\t" + impliedVolLowerQMC + "\t" + impliedVolUpperQMC + "\t");

				/*System.out.println(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t"
						+ pricesMC[m][s] + "\t" + stdErrMC + "\t" + confidenceIntervalMC[0] + "\t" + confidenceIntervalMC[1] + "\t"
						+ pricesQMC[m][s] + "\t" + stdErrQMC + "\t" + confidenceIntervalQMC[0] + "\t" + confidenceIntervalQMC[1] + "\t");*/
			}
		}
	}

	@Test
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
	}
}