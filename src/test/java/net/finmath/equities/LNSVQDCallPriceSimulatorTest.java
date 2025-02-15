package net.finmath.equities;

import net.finmath.equities.models.Black76Model;
import net.finmath.equities.models.LNSVQD.LNSVQDCallPriceSimulator;
import net.finmath.equities.models.LNSVQD.LNSVQDModel;
import net.finmath.equities.models.LNSVQD.LNSVQDModelAnalyticalPricer;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
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
		// Get analytical prices
		double[] pricesAnalytical = lnsvqdModelAnalyticalPricer.getEuropeanOptionPrices(strikeMaturityPairs, true);

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

				List<Integer> seeds = random.ints(5).boxed().collect(Collectors.toList());
				double[] prices = new double[seeds.size()];

				// For statistics
				TDistribution tDistribution = new TDistribution(seeds.size() - 1);

				for(int seed : seeds) {
					LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid);
					lnsvqdCallPriceSimulator.precalculatePaths(seed);
					double simulatedOptionPrice = lnsvqdCallPriceSimulator.getCallPrice(strike, maturity);
					prices[seeds.indexOf(seed)] = simulatedOptionPrice;
				}

				double averagePrice = Arrays.stream(prices).average().getAsDouble();
				pricesMC[m][s] = averagePrice;
				/*double relativeError = Math.abs(averagePrice - lnsvqdOptionValue) / lnsvqdOptionValue;
				double stdError = standardDeviation.evaluate(prices);

				double tQuantile = tDistribution.inverseCumulativeProbability(0.975);
				double lowerConfidenceBound = averagePrice - tQuantile * stdError / seeds.size();
				double upperConfidenceBound = averagePrice + tQuantile * stdError / seeds.size();*/

				/*System.out.println("Average price: " + Arrays.stream(prices).average().getAsDouble() + "; Relative error: " + relativeError
						+ "; Lower 95%-bound: " + lowerConfidenceBound
						+ "; Upper 95%-bound: " + upperConfidenceBound + "\n"
						+ "; Analytical price " + lnsvqdOptionValue
						+ "Analytical price in interval: " + (lowerConfidenceBound <= lnsvqdOptionValue && lnsvqdOptionValue <= upperConfidenceBound));
*/
				/*relativeErrors[m] = relativeError;*/

				System.out.println("ANA: " + pricesAnalytical[m * strikes.length + s] + "\t"
						+ "MC: " + pricesMC[m][s] + "\t"
						+ "BS: " + pricesBS[m][s]);
			}

		}
		/**
		 * Print
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

		System.out.println("BS prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesBS[m][s] + "\t");
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

		LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid);
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