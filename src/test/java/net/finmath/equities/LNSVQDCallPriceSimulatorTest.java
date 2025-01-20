package net.finmath.equities;

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
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.junit.jupiter.api.Assertions.*;

public class LNSVQDCallPriceSimulatorTest {
	/**
	 * Model params
	 */
	// Right params: sigma0=0.8327, theta=1.0139, kappa1=4.8606, kappa2=4.7938, beta=0.1985, volvol=2.3690
	// sigma0=1.5, theta=1.0, kappa1=4.0, kappa2=4.0, beta=0.0, volvol=1.0
	// LogSvParams(sigma0=0.8376, theta=1.0413, kappa1=3.1844, kappa2=3.058, beta=0.1514, volvol=1.8458)
	// LogSvParams(sigma0=0.9778, theta=0.5573, kappa1=4.8360, kappa2=8.6780, beta=2.3128, volvol=1.0484)
	private final double spot0 = 1;
	private final double sigma0 = 0.8376; // 0.8327;
	// Value as in paper
	private final double kappa1 = 3.1844; // 4.8606;
	// Value as in paper
	private final double kappa2 = 3.058; // 4.7938
	private final double theta = 1.0413; // 1.0139
	private final double beta = 0.1514; // 0.1985
	private final double epsilon = 1.8458; // 2.3690;

	/**
	 * Models
	 */
	LNSVQDModel lnsvqdModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);
	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);

	/**
	 * Stat utils
	 */
	StandardDeviation standardDeviation = new StandardDeviation();
	Random random = new Random();

	@Test
	public void getCallPrices() throws CalculationException {
		int numberOfPaths = 100000;
		// Get option values
		double spot = 1;
		double[] maturityGrid = LNSVQDUtils.createTimeGrid(0.2, 1, 4);
		double[] strikes = LNSVQDUtils.createTimeGrid(0.4, 1.4, 2);
		double[] relativeErrors = new double[maturityGrid.length];

		double[][] pricesAnalytical = new double[maturityGrid.length][strikes.length];
		double[][] pricesMC = new double[maturityGrid.length][strikes.length];

		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
					maturity, (int) Math.round(maturity * 365.));
			double riskFreeRate = lnsvqdModelAnalyticalPricer.getRiskFreeRate();
			double discountFactor = Math.exp(-riskFreeRate * maturity);
			double convenienceFcator = 0;
			for(int s = 0; s < strikes.length; s++) {
				double strike = strikes[s];

				double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot, riskFreeRate, sigma0, maturity, strike, true);
				System.out.println("Call oprion price BS: \t" + bsOptionValue);

				double lnsvqdOptionValue = lnsvqdModelAnalyticalPricer.getCallPrice(strike, maturity, discountFactor, convenienceFcator);
				pricesAnalytical[m][s] = lnsvqdOptionValue;
				System.out.println("Call oprion price LNSVQD: \t" + lnsvqdOptionValue);

				List<Integer> seeds = random.ints(5).boxed().collect(Collectors.toList());
				double[] prices = new double[seeds.size()];

				// For statistics
				TDistribution tDistribution = new TDistribution(seeds.size() - 1);

				for(int seed : seeds) {
					LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModel, numberOfPaths, timeGrid);
					// lnsvqdCallPriceSimulator.precalculatePaths(seed);
					double simulatedOptionPrice = 0; // lnsvqdCallPriceSimulator.getCallPrice(strike);
					prices[seeds.indexOf(seed)] = simulatedOptionPrice;
					System.out.println(simulatedOptionPrice);
				}

				double averagePrice = Arrays.stream(prices).average().getAsDouble();
				pricesMC[m][s] = averagePrice;
				double relativeError = Math.abs(averagePrice - lnsvqdOptionValue) / lnsvqdOptionValue;
				double stdError = standardDeviation.evaluate(prices);

				double tQuantile = tDistribution.inverseCumulativeProbability(0.975);
				double lowerConfidenceBound = averagePrice - tQuantile * stdError / seeds.size();
				double upperConfidenceBound = averagePrice + tQuantile * stdError / seeds.size();

				/*System.out.println("Average price: " + Arrays.stream(prices).average().getAsDouble() + "; Relative error: " + relativeError
						+ "; Lower 95%-bound: " + lowerConfidenceBound
						+ "; Upper 95%-bound: " + upperConfidenceBound + "\n"
						+ "; Analytical price " + lnsvqdOptionValue
						+ "Analytical price in interval: " + (lowerConfidenceBound <= lnsvqdOptionValue && lnsvqdOptionValue <= upperConfidenceBound));
*/
				relativeErrors[m] = relativeError;
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
				System.out.print(pricesAnalytical[m][s] + "\t");
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

		/*LNSVQDUtils.printArray(maturityGrid);
		LNSVQDUtils.printArray(relativeErrors);*/
	}

	@Test
	public void printTransformedVolPathMoments() throws CalculationException {
		int numberOfPaths = 100000;
		// Get option values
		double maturity = 1.2;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
				maturity, (int) Math.round(maturity * 365.));

		int seed = 1;

		LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModel, numberOfPaths, timeGrid);
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