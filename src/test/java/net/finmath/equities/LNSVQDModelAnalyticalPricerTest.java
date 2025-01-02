package net.finmath.equities;

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
import org.apache.commons.math3.complex.Complex;
import org.junit.Assert;
import org.junit.jupiter.api.Test;

/**
 * This test
 * 1. compares the semi-analytical LNSVQD call option price to the BS call option price and
 * 2. compares the semi-analytical LNSVQD call option price to the simulated LNSVQD call option price
 */
class LNSVQDModelAnalyticalPricerTest {
	/**
	 * Simulation parameters
	 */
	int numberOfPaths = 1500;
	int seed = 5609;

	/**
	 * Model params
	 */
	private final double spot0 = 1;
	private final double sigma0 = 0.5;
	// Value as in paper
	private final double kappa1 = 2.21;
	// Value as in paper
	private final double kappa2 = 2.18; // 2.18
	private final double theta = 0.4;
	private final double beta = 0.4;
	private final double epsilon = 0.6;

	/**
	 * Models
	 */
	LNSVQDModel lnsvqdSimulationModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);
	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);

	/**
	 * Option params and option
	 */
	double strike = 1;
	double maturity = 1.;
	EuropeanOption europeanOption = new EuropeanOption(maturity, strike, 1, 0);

	/**
	 * Market observables
	 */
	double riskFreeRate = lnsvqdSimulationModel.getRiskFreeRate();
	double discountFactor = Math.exp(-riskFreeRate * maturity);
	double convenienceFcator = 0;

	/**
	 * Time discretization
	 */
	double[] timeGrid = LNSVQDUtils.createTimeGrid(0, 1, 100);
	TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(timeGrid);

	/**
	 * RandomVariableFactory
	 */
	RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

	/**
	 * Tolerance level
	 */
	private final double delta = 10E-3;

	/**
	 * ***************************************************+
	 * 1. Comparison semi-analytical LNSVQD call <-> BS call
	 * ***************************************************+
	 */
	@Test
	void calculateExponentialAffineApproximation() {
		double y = 1.;
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
		Complex exponentialAffineApproximationAnalyticalValue = (charFuncArgs[0].multiply(-lnsvqdModelAnalyticalPricer.getX0())
				.add(charFuncArgs[0].pow(2).add(charFuncArgs[0]).subtract(charFuncArgs[1].multiply(2)).multiply(0.5 * maturity).multiply(Math.pow(lnsvqdModelAnalyticalPricer.getY0(), 2))))
				.exp();
		Complex exponentialAffineApproximationOdeValue = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximation(maturity, charFuncArgs);
		System.out.println("Analytical exponential-affine approximation value at " + maturity + ": " + exponentialAffineApproximationAnalyticalValue);
		System.out.println("ODE-based exponential-affine approximation value at " + maturity + ": " + exponentialAffineApproximationOdeValue);

		Assert.assertEquals(exponentialAffineApproximationAnalyticalValue.getReal(), exponentialAffineApproximationOdeValue.getReal(), delta);
		Assert.assertEquals(exponentialAffineApproximationAnalyticalValue.getImaginary(), exponentialAffineApproximationOdeValue.getImaginary(), delta);
	}

	@Test
	void getCallPrice() {
		// Get option values
		double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot0, riskFreeRate, sigma0, maturity, strike, true);
		double lnsvqdOptionValue = lnsvqdModelAnalyticalPricer.getCallPrice(strike, maturity, discountFactor, convenienceFcator);

		// Print
		System.out.println("Call oprion price BS: \t" + bsOptionValue);
		System.out.println("Call oprion price LNSVQD: \t" + lnsvqdOptionValue);

		Assert.assertEquals(bsOptionValue, lnsvqdOptionValue, delta);
	}

	@Test
	void printODESolution() {
		// Get option values
		double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot0, riskFreeRate, sigma0, maturity, strike, true);
		double lnsvqdOptionValue = lnsvqdModelAnalyticalPricer.getCallPrice(strike, maturity, discountFactor, convenienceFcator);

		// Print
		System.out.println("Call oprion price BS: \t" + bsOptionValue);
		System.out.println("Call oprion price LNSVQD: \t" + lnsvqdOptionValue);

		Assert.assertEquals(bsOptionValue, lnsvqdOptionValue, delta);
	}

	/**
	 * ***************************************************+
	 * 2. Comparison semi-analytical LNSVQD call <-> simulated LNSVQD call
	 * ***************************************************+
	 */

	@Test
	void comparePricingMethods() throws CalculationException {
		// 1. Create the Monte-Carlo Process
		BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
		LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdSimulationModel, brownianMotion);

		/**
		 * NOTE: MonteCarloLNSVQDModel glues together a discretization scheme and a model (e.g. Euler + Heston (terrible idea!));
		 * In our case, the LNSVQD discretization scheme already has an instance of AbstractProcessModel (LNSVQDSimulationModel);
		 * Hence, the model is already "glued" to its discretization scheme and the AssetMonteCarloSimulationModel
		 * is a wrapper; However, we still create a new instance of it in order to be consistent with the finmath-framework
		 */
		MonteCarloLNSVQDModel monteCarloLNSVQDModel = new MonteCarloLNSVQDModel(lnsvqdDiscretizationScheme, seed);

		// Get option values
		double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot0, riskFreeRate, sigma0, maturity, strike, true); // Only valid if volatility is constant!
		double analyticalOptionPrice = lnsvqdModelAnalyticalPricer.getCallPrice(strike, maturity, discountFactor, convenienceFcator);
		double simulatedOptionPrice = europeanOption.getValue(monteCarloLNSVQDModel);

		// Print
		System.out.println("BS option price: \t" + bsOptionValue);
		System.out.println("Semi-analytical option price: \t" + analyticalOptionPrice);
		System.out.println("Simulated option price: \t" + simulatedOptionPrice);

		Assert.assertEquals(analyticalOptionPrice, simulatedOptionPrice, delta);
	}

	/**
	 * ***************************************************+
	 * SECTION 3: Implied SVI surface
	 * ***************************************************+
	 */

	@Test
	public void outputImpliedSVISurface() {
		double endTime = 2;
		int numberOfTimePoints = 10;

		double[] timePoints = LNSVQDUtils.createTimeGrid(0, endTime, numberOfTimePoints - 1);
		double[] moneynessLevels = LNSVQDUtils.createTimeGrid(0.2, 1.8, 8);
		// LNSVQDUtils.printArray(moneynessLevels);
		LNSVQDUtils.printArray(timePoints);

		double[] prices = new double[timePoints.length * moneynessLevels.length];

		for(int k = 0; k < moneynessLevels.length; k++) {
			System.out.print("\t" + moneynessLevels[k]);
		}
		for(int i = 0; i < timePoints.length; i++) {
			System.out.print("\n");
			System.out.print(timePoints[i]);
			for(int j = 0; j < moneynessLevels.length; j++) {
				double discountFactor = Math.exp(-riskFreeRate * timePoints[i]);
				double strike = lnsvqdModelAnalyticalPricer.getSpot0() * moneynessLevels[j];
				double price = lnsvqdModelAnalyticalPricer.getCallPrice(strike, timePoints[i], discountFactor, 0);
				double impliedVol = AnalyticFormulas.blackScholesOptionImpliedVolatility(lnsvqdModelAnalyticalPricer.getSpot0() / discountFactor, timePoints[i], strike, discountFactor, price);
				System.out.print("\t" + impliedVol);
			}
		}
	}

}
