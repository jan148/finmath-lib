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
import org.apache.commons.math3.util.Pair;
import org.junit.Assert;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * This test
 * 1. compares the semi-analytical LNSVQD call option price to the BS call option price and
 * 2. compares the semi-analytical LNSVQD call option price to the simulated LNSVQD call option price
 */
class LNSVQDModelAnalyticalPricerTest extends TestsSetupForLNSVQD{
	/**
	 * Time params
	 */
	LocalDate valuationDate = LocalDate.parse("2000-01-01");

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
	private final double kappa1 = 0; // 3.1844; // 4.8606;
	// Value as in paper
	private final double kappa2 = 0; // 3.058; // 4.7938
	private final double theta = 0; // 1.0413; // 1.0139
	private final double beta = 0; // 0.1514; // 0.1985
	private final double epsilon = 0; // 1.8458; // 2.3690;

	/**
	 * Models
	 */
	LNSVQDModel lnsvqdSimulationModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0, valuationDate);
	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0, valuationDate);

	/**
	 * Option params and option
	 */
	double strike = 1;
	double maturity = 1.2;
	EuropeanOption europeanOption = new EuropeanOption(maturity, strike, 1, 0);

	/**
	 * Market observables
	 */
	// IMPORTANT: Assume constant rf-rate
	double riskFreeRate = lnsvqdSimulationModel.getRiskFreeRate(0);
	double discountFactor = Math.exp(-riskFreeRate * maturity);
	double convenienceFcator = 0;

	/**
	 * Time discretization
	 */
	double[] timeGrid = LNSVQDUtils.createTimeGrid(0, maturity, 100);
	TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(timeGrid);

	/**
	 * RandomVariableFactory
	 */
	RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

	/**
	 * Tolerance level
	 */
	private final double delta = 10E-3;

	@Test
	public void printMatrices() {
		Complex complex0 = Complex.ZERO;
		Complex[][] M0 = {
				{complex0, complex0, complex0, complex0, complex0},
				{complex0, new Complex(0, 0.5), complex0, complex0, complex0},
				{complex0, new Complex(1, 0.5), complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0}};
		System.out.println(M0[2][1]);
	}

	/**
	 * ***************************************************+
	 * 1. Comparison semi-analytical LNSVQD call <-> BS call
	 * ***************************************************+
	 */
	@Test
	public void printAj() {
		int index = 3;
		double ttm = 1;
		double y = 2;
		int numStepsForODEIntegration = (int) (ttm * 365 * lnsvqdModelAnalyticalPricer.numStepsForODEIntegrationPerYear);
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, numStepsForODEIntegration);
		final Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
		Complex[][] solutionPath = lnsvqdModelAnalyticalPricer.getSolutionPathForODESystem(timeGrid, charFuncArgs);
		System.out.println("T \t Real part \t Imaginary part");
		for(int i = 0; i < solutionPath.length; i++) {
			double t = timeGrid[i];
			Complex value = solutionPath[i][index];
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
		}
	}

	@Test
	public void printE2() {
		double ttm = 1;
		double y = 1;
		int numStepsForODEIntegration = (int) (ttm * 365 * lnsvqdModelAnalyticalPricer.numStepsForODEIntegrationPerYear);
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, numStepsForODEIntegration);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximation(t, charFuncArgs).multiply(charFuncArgs[0].multiply(lnsvqdModelAnalyticalPricer.getX0()).add(charFuncArgs[1].multiply(lnsvqdModelAnalyticalPricer.getI0())).exp());
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
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
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
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
	void printODESolution() throws Exception {
		// Get option values
		double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot0, riskFreeRate, sigma0, maturity, strike, true);
		double lnsvqdOptionValue = lnsvqdModelAnalyticalPricer.getCallPrice(strike, maturity);

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
	// IMPORTANT: Assume constant rf-rate
	@Test
	void getCallPrice() throws Exception {
		int numberOfPaths = 50000;
		// Get option values
		double spot = 1;
		double strike = 1.4;
		double maturity = 0.4;
		// IMPORTANT: Assume constant rf-rate
		double discountFactor = Math.exp(-lnsvqdModelAnalyticalPricer.getRiskFreeRate(maturity) * maturity);
		double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot, riskFreeRate, sigma0, maturity, strike, true);
		System.out.println("Call oprion price BS: \t" + bsOptionValue);
		double lnsvqdOptionValue = lnsvqdModelAnalyticalPricer.getCallPrice(strike, maturity);
		System.out.println("Call oprion price LNSVQD: \t" + lnsvqdOptionValue);

		// 1. Create the Monte-Carlo Process
		// List<Integer> seeds = Arrays.asList(1, 2, 3, 4, 5/*, 6, 7, 8, 9, 10*/);
		List<Integer> seeds = Arrays.asList(6, 7, 8, 9, 10);
		// List<Integer> seeds = Arrays.asList(11, 12, 13, 14, 15);
		double[] timeGrid = LNSVQDUtils.createTimeGrid((double) 0,
				maturity, (int) Math.round(maturity * 365.));
		TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(timeGrid);

		/**
		 * NOTE: MonteCarloLNSVQDModel glues together a discretization scheme and a model (e.g. Euler + Heston (terrible idea!));
		 * In our case, the LNSVQD discretization scheme already has an instance of AbstractProcessModel (LNSVQDSimulationModel);
		 * Hence, the model is already "glued" to its discretization scheme and the AssetMonteCarloSimulationModel
		 * is a wrapper; However, we still create a new instance of it in order to be consistent with the finmath-framework
		 */
		// Get option values
		// double[] values = new double[seeds.length];
		System.out.println("Available processors: " + Runtime.getRuntime().availableProcessors());
		ForkJoinPool forkJoinPool = new ForkJoinPool(4);

		Function<Integer, Double> getCallValue = seed -> {
			String result = "Seed " + seed + " processed by " + Thread.currentThread().getName();
			try {
				// System.out.println(result);
				BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
				LNSVQDModel lnsvqdSimulationModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0, valuationDate);
				LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdSimulationModel, brownianMotion);
				MonteCarloLNSVQDModel monteCarloLNSVQDModel = new MonteCarloLNSVQDModel(lnsvqdDiscretizationScheme, seed);
				EuropeanOption europeanOption = new EuropeanOption(maturity, strike, 1, 0);
				double simulatedOptionPrice = europeanOption.getValue(monteCarloLNSVQDModel);
				System.out.println("Seed " + seed + ": \t" + simulatedOptionPrice);
			} catch(CalculationException e) {
				Thread.currentThread().interrupt();
			}
			return 0.;
		};

		List<Double> results = seeds.parallelStream()
				.map(getCallValue) // Apply the function to each seed
				.collect(Collectors.toList());
	}

	/*@Test
	public void getEuropeanOptionPrices() throws Exception {
		double spot = 1;
		double[] maturityGrid = LNSVQDUtils.createTimeGrid(0.2, 1, 4);
		double[] strikes = LNSVQDUtils.createTimeGrid(0.6, 1.4, 4);

		*//**
		 * Get all the integration points from the integrator, in our case Simpson
		 *//*
		List<Double> yGridForIntegration = new ArrayList<>();

		// Next lines adapted from finmath's Simpson implementation
		final double	lowerBound			= lnsvqdModelAnalyticalPricer.lowerBound;
		final double	upperBound			= lnsvqdModelAnalyticalPricer.upperBound;
		final double	range				= upperBound-lowerBound;

		final int  numberOfEvaluationPoints = (int) upperBound * 10; // Need to change this in accordance with LNSVQD pricer
		final int numberOfDoubleSizeIntervals	= (int) ((numberOfEvaluationPoints-1) / 2.0);

		final double doubleInterval = range / numberOfDoubleSizeIntervals;
		final double singleInterval = 0.5 * doubleInterval;

		IntStream intervals = IntStream.range(1, numberOfDoubleSizeIntervals);

		intervals.forEach(
				i -> {
					yGridForIntegration.add(lowerBound + i * doubleInterval);
					yGridForIntegration.add(lowerBound + i * doubleInterval + singleInterval);
				}
		);

		yGridForIntegration.add(lowerBound + singleInterval);
		yGridForIntegration.add(lowerBound);
		yGridForIntegration.add(upperBound);

		List<Double> yGridForIntegrationSortedDistinct = yGridForIntegration.stream().sorted().distinct().collect(Collectors.toList());
		*//**
		 * Get call prices
		 *//*
		double[] callPricesMethod1 = lnsvqdModelAnalyticalPricer.getEuropeanOptionPrices(strikes, maturityGrid, yGridForIntegrationSortedDistinct);
		List<Pair<Double, Double>> mesh = LNSVQDUtils.create2dMesh(maturityGrid, strikes);
		double[] callPricesMethod2 = lnsvqdModelAnalyticalPricer.getEuropeanOptionPrices(mesh);

		*//**
		 * Print
		 *//*
		LNSVQDUtils.printPricesFromMaturityStrikeGrid(maturityGrid, strikes, callPricesMethod1);
		LNSVQDUtils.printPricesFromMaturityStrikeGrid(maturityGrid, strikes, callPricesMethod2);

	}*/

	/**
	 * ***************************************************+
	 * SECTION 3: TEST CALL PRICE CALCULATION WITH PRECALCULATION OF E2
	 * ***************************************************+
	 */
	@Test
	public void calculateExponentialAffineApproximationFullPathTest() {
		double ttm = 1;
		double y = 2.;
		int numStepsForODEIntegration = (int) (ttm * 365 * lnsvqdModelAnalyticalPricer.numStepsForODEIntegrationPerYear);
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, numStepsForODEIntegration);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		Complex[] values = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximationFullPath(timeGrid, charFuncArgs);
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = values[i].multiply(charFuncArgs[0].multiply(lnsvqdModelAnalyticalPricer.getX0()).add(charFuncArgs[1].multiply(lnsvqdModelAnalyticalPricer.getI0())).exp());
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
		}
	}

	/**
	 * ***************************************************+
	 * SECTION 4: Implied SVI surface
	 * ***************************************************+
	 */
	@Test
	public void outputImpliedSVISurface() throws Exception {
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
				double price = lnsvqdModelAnalyticalPricer.getCallPrice(strike, timePoints[i]);
				double impliedVol = AnalyticFormulas.blackScholesOptionImpliedVolatility(lnsvqdModelAnalyticalPricer.getSpot0() / discountFactor, timePoints[i], strike, discountFactor, price);
				System.out.print("\t" + impliedVol);
			}
		}
	}

	/**
	 * ***************************************************+
	 * SECTION 4: TEST CALL PRICE CALCULATION WITH PRECALCULATION OF E2
	 * ***************************************************+
	 */

	@Test
	public void printPricesWrtPricerParams() throws Exception {
		double a;
		ArrayList<Pair<Double, Double>> strikeMatPairs = setDAXHestonSetupSIM();

		double[] volAna = lnsvqdModelAnalyticalPricer.getImpliedVolsStrikeMatList(strikeMatPairs, null);

	}
}
