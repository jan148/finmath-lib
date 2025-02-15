package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.LNSVQDCallPriceSimulator;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
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

		// Get analytical prices
		double[] pricesAnalytical = lnsvqdModelAnalyticalPricer.getEuropeanOptionPrices(strikeMaturityPairs, true);

		// Get MC-prices
		List<Integer> seeds = random.ints(5).boxed().collect(Collectors.toList());
		// seeds * maturities * strikes
		double[][][] pricesMCPerSeed = new double[seeds.size()][maturityGrid.length][strikes.length];

		for(int j = 0; j < seeds.size(); j++) {
			int seed = seeds.get(j);
			lnsvqdCallPriceSimulator.precalculatePaths(seed);
			for(int m = 0; m < maturityGrid.length; m++) {
				double maturity = maturityGrid[m];
				for(int s = 0; s < strikes.length; s++) {
					double strike = strikes[s];
					pricesMCPerSeed[j][m][s] = lnsvqdCallPriceSimulator.getCallPrice(strike, maturity);
				}
			}
		}

		double[][] pricesMC = new double[maturityGrid.length][strikes.length];
		for(int m = 0; m < maturityGrid.length; m++) {
			int finalM = m;
			for(int s = 0; s < strikes.length; s++) {
				int finalS = s;
				pricesMC[m][s] = IntStream.range(0, seeds.size())
						.mapToDouble(i -> pricesMCPerSeed[i][finalM][finalS])
						.average()
						.getAsDouble();
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