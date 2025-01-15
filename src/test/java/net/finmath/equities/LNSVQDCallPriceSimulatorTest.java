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
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.*;

public class LNSVQDCallPriceSimulatorTest {
	/**
	 * Model params
	 */
	// Right params: sigma0=0.8327, theta=1.0139, kappa1=4.8606, kappa2=4.7938, beta=0.1985, volvol=2.3690
	private final double spot0 = 1;
	private final double sigma0 = 0.8327; //0.41;
	// Value as in paper
	private final double kappa1 =  4.8606; // 4.8606;
	// Value as in paper
	private final double kappa2 = 4.7938; // 4.7938
	private final double theta =  1.0139; // 1.0139
	private final double beta = 0.1985; // 0.1985
	private final double epsilon = 2.3690; // 2.3690;

	/**
	 * Models
	 */
	LNSVQDModel lnsvqdModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);
	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);

	@Test
	public void getCallPrice() throws CalculationException {
		int numberOfPaths = 30000;
		// Get option values
		double spot = 1;
		double strike = 1.4;
		double maturity = 0.4;
		double[] timeGrid = LNSVQDUtils.createTimeGrid((double) 0,
				maturity, (int) Math.round(maturity * 365.));
		double riskFreeRate = lnsvqdModelAnalyticalPricer.getRiskFreeRate();
		double discountFactor = Math.exp(-riskFreeRate * maturity);

		double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot, riskFreeRate, sigma0, maturity, strike, true);
		System.out.println("Call oprion price BS: \t" + bsOptionValue);

		double convenienceFcator = 0;
		double lnsvqdOptionValue = lnsvqdModelAnalyticalPricer.getCallPrice(strike, maturity, discountFactor, convenienceFcator);
		System.out.println("Call oprion price LNSVQD: \t" + lnsvqdOptionValue);

		List<Integer> seeds = Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
		double[] prices = new double[seeds.size()];

		for(int seed : seeds) {
			LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModel, numberOfPaths, timeGrid);
			lnsvqdCallPriceSimulator.precalculatePaths(seed);
			double simulatedOptionPrice = lnsvqdCallPriceSimulator.getCallPrice(strike);
			prices[seeds.indexOf(seed)] = simulatedOptionPrice;
			System.out.println(seed + ": \t" + simulatedOptionPrice);
		}

		System.out.println("Average price: " + Arrays.stream(prices).average().getAsDouble());

	}
}