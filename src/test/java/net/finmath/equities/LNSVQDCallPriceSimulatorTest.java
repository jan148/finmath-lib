package net.finmath.equities;

import net.finmath.equities.models.Black76Model;
import net.finmath.equities.models.LNSVQD.*;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.*;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloLNSVQDModel;
import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.montecarlo.process.LNSVQDDiscretizationScheme;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import java.time.LocalDate;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.junit.jupiter.api.Assertions.*;

public class LNSVQDCallPriceSimulatorTest extends TestsSetupForLNSVQD {
	/**
	 * Stat utils
	 */
	StandardDeviation standardDeviation = new StandardDeviation();
	Random random = new Random();

	@Test
	public void getCallPrices() throws Exception {
		// Get option values
		double[] strikes = LNSVQDUtils.createTimeGrid(0.6 * spot0, 1.4 * spot0, 4);
		List<Pair<Double, Double>> strikeMaturityPairs = LNSVQDUtils.create2dMesh(maturityGrid, strikes);
		//double[] relativeErrors = new double[maturityGrid.length];

		double[][] pricesBS = new double[maturityGrid.length][strikes.length];
		double[][] pricesMC = new double[maturityGrid.length][strikes.length];
		double[][] pricesQMC = new double[maturityGrid.length][strikes.length];
		double[][] pricesMCFinmath = new double[maturityGrid.length][strikes.length];
		// Get analytical prices
		double[] pricesAnalytical = lnsvqdModelAnalyticalPricer.getEuropeanOptionPrices(strikeMaturityPairs, true);

		double[][] stdErrorsMc = new double[maturityGrid.length][strikes.length];
		double[][] stdErrorsMcFinmath = new double[maturityGrid.length][strikes.length];
		double[][] stdErrorsQMc = new double[maturityGrid.length][strikes.length];

		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			double discountFactor = equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
			double forward = spot0 / discountFactor;
			double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
					maturity, (int) Math.round(maturity * 365.));
			for(int s = 0; s < strikes.length; s++) {
				double strike = strikes[s];

				// BS value
				pricesBS[m][s] = Black76Model.optionPrice(forward, strike, maturity, selectedParams[0], true, discountFactor);

				List<Integer> seeds = random.ints(10).boxed().collect(Collectors.toList());
				double[] prices = new double[seeds.size()];
				double[] pricesFm = new double[seeds.size()];
				double[] pricesQ = new double[seeds.size()];

				for(int seed : seeds) {
					// Normal MC
					LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdCallPriceSimulator.precalculatePaths(seed);
					double simulatedOptionPrice = lnsvqdCallPriceSimulator.getCallPrice(strike, maturity);
					prices[seeds.indexOf(seed)] = simulatedOptionPrice;

					// QMC
					LNSVQDPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
					double simulatedOptionPriceQMC = lnsvqdPriceSimulatorQMC.getCallPrice(strike, maturity);
					pricesQ[seeds.indexOf(seed)] = simulatedOptionPriceQMC;

					// Price for finmath implementation
					// BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
					BrownianMotionFromMersenneRandomNumbers brownianBridge = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
					LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdModelAnalyticalPricer, brownianBridge);
					MonteCarloLNSVQDModel monteCarloLNSVQDModel = new MonteCarloLNSVQDModel(lnsvqdDiscretizationScheme, seed);
					EuropeanOption europeanOption = new EuropeanOption(maturity, strike, 1, 0);
					pricesFm[seeds.indexOf(seed)] = 123456789; // europeanOption.getValue(monteCarloLNSVQDModel);
				}

				double averagePrice = Arrays.stream(prices).average().getAsDouble();
				double varMC = Arrays.stream(prices).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
				double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());

				double averagePriceFm = Arrays.stream(pricesFm).average().getAsDouble();
				double varMcFm = Arrays.stream(pricesFm).map(x -> Math.pow(x - averagePriceFm, 2)).sum() / (seeds.size() - 1);
				double stdErrMcFm = Math.sqrt(varMcFm) / Math.sqrt(seeds.size());

				double averagePriceQMC = Arrays.stream(pricesQ).average().getAsDouble();
				double varQMC = Arrays.stream(pricesQ).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
				double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());

				pricesMC[m][s] = averagePrice;
				pricesMCFinmath[m][s] = averagePriceFm;
				pricesQMC[m][s] = averagePriceQMC;

				stdErrorsMc[m][s] = stdErrMC;
				stdErrorsMcFinmath[m][s] = stdErrMcFm;
				stdErrorsQMc[m][s] = stdErrQMC;

				System.out.println("ANA: " + pricesAnalytical[m * strikes.length + s] + "\t"
						+ "MC: " + pricesMC[m][s] + " (" + stdErrMC + ")" + "\t"
						+ "MC finmath: " + pricesMCFinmath[m][s] + " (" + stdErrMcFm + ")" + "\t"
						+ "QMC: " + pricesQMC[m][s] + " (" + stdErrQMC + ")" +"\t"
						+ "BS: " + pricesBS[m][s]);
			}

		}
		/**
		 * Print prices
		 */
		System.out.println("Analytical prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesAnalytical[m * strikes.length + s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("MC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesQMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("BS prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesBS[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("Mc prices finmath");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesMCFinmath[m][s] + "\t");
			}
			System.out.print("\n");
		}

		/**
		 * Print standard error
		 */
		System.out.println("MC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(stdErrorsMc[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(stdErrorsQMc[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("MC finmath std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(stdErrorsMcFinmath[m][s] + "\t");
			}
			System.out.print("\n");
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

		LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
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