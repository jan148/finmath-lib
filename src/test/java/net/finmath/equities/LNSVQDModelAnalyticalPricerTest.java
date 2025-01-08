package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.ComplexRungeKutta4thOrderIntegrator;
import net.finmath.equities.models.LNSVQD.LNSVQDModel;
import net.finmath.equities.models.LNSVQD.LNSVQDModelAnalyticalPricer;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.integration.SimpsonRealIntegrator;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloLNSVQDModel;
import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.montecarlo.process.LNSVQDDiscretizationScheme;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.complex.Complex;
import org.junit.Assert;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.DoubleUnaryOperator;

/**
 * This test
 * 1. compares the semi-analytical LNSVQD call option price to the BS call option price and
 * 2. compares the semi-analytical LNSVQD call option price to the simulated LNSVQD call option price
 */
class LNSVQDModelAnalyticalPricerTest {
	/**
	 * Simulation parameters
	 */
	int numberOfPaths = 20000;
	int seed = 3;

	/**
	 * Model params
	 */
	// Right params: sigma0=0.8327, theta=1.0139, kappa1=4.8606, kappa2=4.7938, beta=0.1985, volvol=2.3690
	private final double spot0 = 1;
	private final double sigma0 = 0.8327; //0.41;
	// Value as in paper
	private final double kappa1 = 4.8606; //2.21;
	// Value as in paper
	private final double kappa2 = 4.7938; //2.18;
	private final double theta = 1.0139; //0.38;
	private final double beta = 0.1985; //0.5;
	private final double epsilon = 2.3690; //3.06;

	/**
	 * Models
	 */
	LNSVQDModel lnsvqdSimulationModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);
	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);

	/**
	 * Option params and option
	 */
	double strike = 1;
	double maturity = 1 / 60.;
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
	double[] timeGrid = LNSVQDUtils.createTimeGrid(0, maturity, 99);
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
	public void delete() {
		System.out.println(4.8606 * 1.0139);
		System.out.println(2.21 * 0.38);
	}
	@Test
	public void printAj() {
		int index = 4;
		double ttm = 0;
		double y = 30;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, lnsvqdModelAnalyticalPricer.numStepsForODEIntegration);
		final Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
		Complex[][] solutionPath = lnsvqdModelAnalyticalPricer.getSolutionPathForODESystem(timeGrid, charFuncArgs);
		System.out.println("T \t Real part \t Imaginary part");
		for(int i = 0; i < solutionPath.length; i++) {
			double t = timeGrid[i];
			Complex value = solutionPath[i][index];
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t"+  imaginaryPart);
		}
	}

	@Test
	public void printE2() {
		double ttm = 1;
		double y = 200;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, lnsvqdModelAnalyticalPricer.numStepsForODEIntegration);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximation(t, charFuncArgs).multiply(charFuncArgs[0].multiply(lnsvqdModelAnalyticalPricer.getX0()).add(charFuncArgs[1].multiply(lnsvqdModelAnalyticalPricer.getI0())).exp());
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t"+  imaginaryPart);
		}
	}

	@Test
	public void printMGF() {
		double ttm = 0;
		double y = 2;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, 100);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximation(t, charFuncArgs);
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t"+  imaginaryPart);
		}
	}

	// Next test work only if volatility parameters are 0, i.e. model is a GBM
	@Test
	void calculateExponentialAffineApproximation() {
		if(kappa1 != 0 || kappa2 != 0 || theta != 0 || beta != 0 || epsilon != 0) {
			return;
		}
		double y = 1.;
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
		Complex exponentialAffineApproximationAnalyticalValue = (charFuncArgs[0].multiply(-lnsvqdModelAnalyticalPricer.getX0())
				.add(charFuncArgs[0].pow(2).add(charFuncArgs[0]).subtract(charFuncArgs[1].multiply(2)).multiply(0.5 * maturity)
						.multiply(Math.pow(lnsvqdModelAnalyticalPricer.getY0(), 2)))).exp();
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
	public void testPriceAt0() {
		double ttm = 0;
		double strike = 1;
		double logMoneyness = Math.log(lnsvqdModelAnalyticalPricer.getSpot0() / strike);
		double discountFactor = 1;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, 100);

		/*DoubleUnaryOperator integrand = new DoubleUnaryOperator() {
			@Override
			public double applyAsDouble(double y) {
				// 1. Compute the value of the affine-exponential approximation
				double E2 = 1;

				// 2. Calculate result
				double result = (1 / (y * y + 0.25)) * E2;
				return result;
			}
		};*/

		UnivariateFunction function = x -> (1 / (x * x + 0.25));

		// Create an integrator instance
		// Parameters: number of points, relative accuracy, absolute accuracy, min and max iterations
		int numberOfPoints = 10; // Number of integration points for Legendre-Gauss quadrature
		double relativeAccuracy = 1.0e-6;
		double absoluteAccuracy = 1.0e-9;
		int minIterations = 3;
		int maxIterations = 100;

		IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
				numberOfPoints, relativeAccuracy, absoluteAccuracy, minIterations, maxIterations
		);

		// Define the integration interval [a, b]
		double lowerBound = 0;
		double upperBound = 1000000;

		// Perform the integration
		double result = 0;
		try {
			result = integrator.integrate(100000000, function, lowerBound, upperBound);
			System.out.println("Integral result: " + result);
		} catch (Exception e) {
			System.err.println("Integration failed: " + e.getMessage());
		}
		double optionPrice = lnsvqdModelAnalyticalPricer.getSpot0() - (discountFactor * strike / Math.PI) * result;
		System.out.println(optionPrice);

		/**
		 * Integerate
		 */
		/*Complex[] state = new Complex[]{new Complex(0., 0.)};
		List<BiFunction<Double, Complex[], Complex>> odeSystem = new ArrayList<>();
		odeSystem.add(integrand);
		ComplexRungeKutta4thOrderIntegrator complexRungeKutta4thOrderIntegrator = new ComplexRungeKutta4thOrderIntegrator(state, odeSystem);
		// Choose the end point of the solution path
		Complex integral = complexRungeKutta4thOrderIntegrator.getSolutionPath
				(lnsvqdModelAnalyticalPricer.yGridForInfiniteIntegral)[lnsvqdModelAnalyticalPricer.numStepsForInfiniteIntegral][0];

		double integralReal = integral.getReal();
		double optionPrice = lnsvqdModelAnalyticalPricer.getSpot0() - (discountFactor * strike / Math.PI) * integralReal;

		System.out.println(optionPrice);*/
	}

	@Test
	void simulationPrice() throws CalculationException {
		long startTime = System.nanoTime();

		{
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
			double simulatedOptionPrice = europeanOption.getValue(monteCarloLNSVQDModel);

			// Print
			System.out.println("Simulated option price: \t" + simulatedOptionPrice);
		}

		long endTime = System.nanoTime(); // Record end time
		long duration = endTime - startTime; // Calculate duration in nanoseconds

		System.out.println("Execution time: " + duration / 1_000_000 + " milliseconds.");
	}

	@Test
	void comparePricingMethods() throws CalculationException {
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0.0, 1.2, 5); // Change back to 0.0 = t0
		double[] strikes = LNSVQDUtils.createTimeGrid(0.4, 1.6, 2);
		for (int i = 0; i < strikes.length; i++) {
			strikes[i] *= spot0;
		}

		// 1. Create the Monte-Carlo Process
		TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(timeGrid);
		BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
		LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdSimulationModel, brownianMotion);

		/**
		 * NOTE: MonteCarloLNSVQDModel glues together a discretization scheme and a model (e.g. Euler + Heston (terrible idea!));
		 * In our case, the LNSVQD discretization scheme already has an instance of AbstractProcessModel (LNSVQDSimulationModel);
		 * Hence, the model is already "glued" to its discretization scheme and the AssetMonteCarloSimulationModel
		 * is a wrapper; However, we still create a new instance of it in order to be consistent with the finmath-framework
		 */
		MonteCarloLNSVQDModel monteCarloLNSVQDModel = new MonteCarloLNSVQDModel(lnsvqdDiscretizationScheme, seed);

		// 1. Print analytical prices
		/*for(int j = 0; j < strikes.length; j++) {
			System.out.print("\t" + strikes[j]);
		}
		System.out.print("\n");
		for(int i = 0; i < timeGrid.length; i++) {
			double ttm = timeGrid[i];
			System.out.print(ttm);
			for(int j = 0; j < strikes.length; j++) {
				double strike = strikes[j];
				// Get option values
				double analyticalOptionPrice = lnsvqdModelAnalyticalPricer.getCallPrice(strike, ttm, discountFactor, convenienceFcator);
				System.out.print("\t" + analyticalOptionPrice);
			}
			System.out.print("\n");
		}*/

		// 2. Print simulation prices
		for(int j = 0; j < strikes.length; j++) {
			System.out.print("\t" + strikes[j]);
		}
		System.out.print("\n");
		for(int i = 0; i < timeGrid.length; i++) {
			double ttm = timeGrid[i];
			System.out.print(ttm);
			for(int j = 0; j < strikes.length; j++) {
				double strike = strikes[j];
				// Get option values
				EuropeanOption europeanOption = new EuropeanOption(ttm, strike, 1, 0);
				double simulatedOptionPrice = europeanOption.getValue(0., monteCarloLNSVQDModel).getAverage();
				System.out.print("\t" + simulatedOptionPrice);
			}
			System.out.print("\n");
		}
	}

	/**
	 * ***************************************************+
	 * SECTION 3: TEST CALL PRICE CALCULATION WITH PRECALCULATION OF E2
	 * ***************************************************+
	 */
	@Test
	public void calculateExponentialAffineApproximationFullPathTest() {
		double ttm = 1;
		double y = 2.;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, lnsvqdModelAnalyticalPricer.numStepsForODEIntegration);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		Complex[] values = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximationFullPath(timeGrid, charFuncArgs);
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = values[i].multiply(charFuncArgs[0].multiply(lnsvqdModelAnalyticalPricer.getX0()).add(charFuncArgs[1].multiply(lnsvqdModelAnalyticalPricer.getI0())).exp());
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t"+  imaginaryPart);
		}
	}

	@Test
	public void comparePricingFunctions() {
		double[] timeGrid = {0.2};
		double[] strikes = {lnsvqdModelAnalyticalPricer.getSpot0()};

		double optionPrice = lnsvqdModelAnalyticalPricer.getCallPrice(strikes[0], timeGrid[0], discountFactor, 0);
		double[] optionPrices;
		try {
			optionPrices = lnsvqdModelAnalyticalPricer.getCallPrices(strikes, timeGrid);
		} catch(Exception e) {
			throw new IllegalArgumentException(e);
		}

		System.out.println("Price old method: " + optionPrice + "\tPrice new method: " + optionPrices[0]);
	}

	@Test
	public void getCallPricesTest() {
		double[] timeGrid = {0}; //LNSVQDUtils.createTimeGrid(0.0, 1.2, 5); // Change back to 0.0 = t0
		double[] strikes = LNSVQDUtils.createTimeGrid(0.4, 1.6, 2);
		for (int i = 0; i < strikes.length; i++) {
			strikes[i] *= spot0;
		}

		double[] optionPrices;
		try {
			optionPrices = lnsvqdModelAnalyticalPricer.getCallPrices(strikes, timeGrid);
		} catch(Exception e) {
			throw new IllegalArgumentException(e);
		}

		for(int k = 0; k < strikes.length; k++) {
			System.out.print("\t" + strikes[k]);
		}
		for(int i = 0; i < timeGrid.length; i++) {
			System.out.print("\n");
			System.out.print(timeGrid[i]);
			for(int j = 0; j < strikes.length; j++) {
				double discountFactor = Math.exp(-riskFreeRate * timeGrid[i]);
				double strike = strikes[j];
				double price = optionPrices[i * strikes.length + j];
				//double impliedVol = AnalyticFormulas.blackScholesOptionImpliedVolatility(lnsvqdModelAnalyticalPricer.getSpot0() / discountFactor, timePoints[i], strike, discountFactor, price);
				System.out.print("\t" + price);
			}
		}
	}

	/**
	 * ***************************************************+
	 * SECTION 4: Implied SVI surface
	 * ***************************************************+
	 */
	@Test
	public void outputImpliedSVISurface() {
		double endTime = 1;
		int numberOfTimePoints = 10;

		double[] timePoints = LNSVQDUtils.createTimeGrid(0, endTime, numberOfTimePoints - 1);
		double[] moneynessLevels = LNSVQDUtils.createTimeGrid(0.6, 1.4, 8);
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
