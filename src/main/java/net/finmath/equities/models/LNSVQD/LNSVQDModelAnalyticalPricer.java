package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.Black76Model;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.integration.SimpsonRealIntegrator;
import net.finmath.time.daycount.DayCountConvention;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Pair;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class LNSVQDModelAnalyticalPricer extends LNSVQDModel {
	/**
	 * Numerical parameters
	 */
	// 1. For ODE-solution
	public int numStepsForODEIntegrationPerYear = 365;

	// 2. For unbounded integration
	// Integration bounds params
	public final double lowerBound = 0;
	public double upperBound = 100; // 100;
	List<Double> yGridForIntegration = new ArrayList<>();

	// finmath integrator
	int numberOfEvaluationPoints = (int) upperBound * 10;
	SimpsonRealIntegrator simpsonRealIntegrator = new SimpsonRealIntegrator(lowerBound, upperBound, numberOfEvaluationPoints, false);

	public LNSVQDModelAnalyticalPricer(double spot0, double sigma0, double kappa1, double kappa2, double theta, double beta, double epsilon, double I0, LocalDate valuationDate, EquityForwardStructure equityForwardStructure) {
		super(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0, valuationDate, equityForwardStructure);

		/**
		 * Get all the integration points from the integrator, in our case Simpson
		 */
		// Next lines adapted from finmath's Simpson implementation
		final double range = upperBound - lowerBound;

		final int numberOfDoubleSizeIntervals = (int) ((numberOfEvaluationPoints - 1) / 2.0);

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

		yGridForIntegration = yGridForIntegration.stream().sorted().distinct().collect(Collectors.toList());
	}

	public LNSVQDModelAnalyticalPricer(double spot0, double sigma0, double kappa1, double kappa2, double theta, double beta, double epsilon, double I0, LocalDate valuationDate) {
		super(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0, valuationDate);

		/**
		 * Get all the integration points from the integrator, in our case Simpson
		 */
		// Next lines adapted from finmath's Simpson implementation
		final double range = upperBound - lowerBound;

		final int numberOfDoubleSizeIntervals = (int) ((numberOfEvaluationPoints - 1) / 2.0);

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

		yGridForIntegration = yGridForIntegration.stream().sorted().distinct().collect(Collectors.toList());
	}

	public void setUpperBoundForIntegration(double upperBoundForIntegration) {
		numberOfEvaluationPoints = (int) upperBoundForIntegration * 10;
		simpsonRealIntegrator = new SimpsonRealIntegrator(lowerBound, upperBoundForIntegration, numberOfEvaluationPoints, false);

		/**
		 * Get all the integration points from the integrator, in our case Simpson
		 */
		// Next lines adapted from finmath's Simpson implementation
		final double range = upperBoundForIntegration - lowerBound;

		final int numberOfDoubleSizeIntervals = (int) ((numberOfEvaluationPoints - 1) / 2.0);

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
		yGridForIntegration.add(upperBoundForIntegration);

		yGridForIntegration = yGridForIntegration.stream().sorted().distinct().collect(Collectors.toList());
	}

	/**
	 * ***************************************************+
	 * SECTION 1: Semi-analytical call option price calculation
	 * ***************************************************+
	 */
	public ArrayList<BiFunction<Double, Complex[], Complex>> getOdeSystem(Complex[] charFuncArgs) {
		/**
		 * Create the functions for the exponential-affine approximation; Names like in the paper
		 */
		// Some constants
		final Complex complex0 = Complex.ZERO;
		final Complex totalInstVar = new Complex(this.totalInstVar, 0);
		final Complex mixedDeg3 = new Complex(this.theta * this.totalInstVar, 0);
		final Complex mixedDeg4 = new Complex(Math.pow(this.theta, 2) * this.totalInstVar, 0);

		// 1. Matrices
		Complex[][] M0 = { { complex0, complex0, complex0, complex0, complex0 },
				{ complex0, mixedDeg4.multiply(0.5), complex0, complex0, complex0 },
				{ complex0, complex0, complex0, complex0, complex0 },
				{ complex0, complex0, complex0, complex0, complex0 },
				{ complex0, complex0, complex0, complex0, complex0 } };

		Complex[][] M1 = { { complex0, complex0, complex0, complex0, complex0 },
				{ complex0, mixedDeg3, mixedDeg4, complex0, complex0 },
				{ complex0, mixedDeg4, complex0, complex0, complex0 },
				{ complex0, complex0, complex0, complex0, complex0 },
				{ complex0, complex0, complex0, complex0, complex0 } };

		Complex[][] M2 = { { complex0, complex0, complex0, complex0, complex0 },
				{ complex0, totalInstVar.multiply(1. / 2), mixedDeg3.multiply(2), mixedDeg4.multiply(3. / 2),
						complex0 },
				{ complex0, mixedDeg3.multiply(2), mixedDeg4.multiply(2), complex0, complex0 },
				{ complex0, mixedDeg4.multiply(3. / 2), complex0, complex0, complex0 },
				{ complex0, complex0, complex0, complex0, complex0 } };

		Complex[][] M3 = { { complex0, complex0, complex0, complex0, complex0 },
				{ complex0, complex0, totalInstVar, mixedDeg3.multiply(3), mixedDeg4.multiply(2) },
				{ complex0, totalInstVar, mixedDeg3.multiply(4), mixedDeg4.multiply(3), complex0 },
				{ complex0, mixedDeg3.multiply(3), mixedDeg4.multiply(3), complex0, complex0 },
				{ complex0, mixedDeg4.multiply(2), complex0, complex0, complex0 } };

		Complex[][] M4 = { { complex0, complex0, complex0, complex0, complex0 },
				{ complex0, complex0, complex0, totalInstVar.multiply(3. / 2), mixedDeg3.multiply(4) },
				{ complex0, complex0, totalInstVar.multiply(2), mixedDeg3.multiply(6), mixedDeg4.multiply(4) },
				{ complex0, totalInstVar.multiply(3. / 2), mixedDeg3.multiply(6), mixedDeg4.multiply(9. / 2),
						complex0 },
				{ complex0, mixedDeg3.multiply(4), mixedDeg4.multiply(4), complex0, complex0 } };

        /*# fills Ls
        L = np.zeros((n, n), dtype=np.complex128)
        L[0, 1], L[0, 2] = lamda - theta2 * beta * vol_backbone_eta * phi, qv2
        L[1, 1], L[1, 2] = -kappa_p - 2.0 * theta * beta * vol_backbone_eta * phi, 2.0 * (lamda + qv - theta2 * beta * vol_backbone_eta * phi)
        L[2, 1], L[2, 2] = -kappa2_p - beta * vol_backbone_eta * phi, vartheta2 - 2.0 * kappa_p - 4.0 * theta * beta * vol_backbone_eta * phi

        if expansion_order == ExpansionOrder.SECOND:
        L[1, 3] = 3.0*qv2
        L[2, 3], L[2, 4] = 3.0 * (2.0 * qv - theta2 * beta * vol_backbone_eta * phi), 6.0 * qv2
        L[3, 2], L[3, 3], L[3, 4] = -2.0 * (kappa2_p + beta * vol_backbone_eta * phi), 3.0 * (vartheta2 - kappa_p - 2.0 * theta * beta * vol_backbone_eta * phi), 4.0 * (3.0 * qv - theta2 * beta * vol_backbone_eta * phi)
        L[4, 3], L[4, 4] = -3.0 * (kappa2_p + beta * vol_backbone_eta * phi), 2.0 * (vartheta2 - 2.0 * kappa_p - 4.0 * theta * beta * vol_backbone_eta * phi)*/

		// 2. Vectors
		Complex L01 = complex0;
		Complex L02 = charFuncArgs[0].multiply(-Math.pow(theta, 2) * beta);
		Complex L03 = mixedDeg4;
		Complex L04 = complex0;
		Complex L05 = complex0;
		Complex[] L0 = { L01, L02, L03, L04, L05 };

		Complex L11 = complex0;
		Complex L12 = charFuncArgs[0].multiply(-theta * beta * 2).subtract(kappa1 + kappa2 * theta);
		Complex L13 = mixedDeg3.subtract(charFuncArgs[0].multiply(Math.pow(theta, 2) * beta)).multiply(2);
		Complex L14 = mixedDeg4.multiply(3);
		Complex L15 = complex0;
		Complex[] L1 = { L11, L12, L13, L14, L15 };
		// L[2, 3], L[2, 4] = 3.0 * (2.0 * qv - theta2 * beta * vol_backbone_eta * phi), 6.0 * qv2
		Complex L21 = complex0;
		Complex L22 = charFuncArgs[0].multiply(-beta).subtract(kappa2);
		Complex L23 = totalInstVar.subtract(2 * (kappa1 + kappa2 * theta))
				.subtract(charFuncArgs[0].multiply(4 * theta * beta));
		Complex L24 = mixedDeg3.multiply(2).subtract(charFuncArgs[0].multiply(Math.pow(theta, 2) * beta)).multiply(3); // !
		Complex L25 = mixedDeg4.multiply(6);
		Complex[] L2 = { L21, L22, L23, L24, L25 };
		// L[3, 2], L[3, 3], L[3, 4] = -2.0 * (kappa2_p + beta * vol_backbone_eta * phi), 3.0 * (vartheta2 - kappa_p - 2.0 * theta * beta * vol_backbone_eta * phi), 4.0 * (3.0 * qv - theta2 * beta * vol_backbone_eta * phi)
		Complex L31 = complex0;
		Complex L32 = complex0;
		Complex L33 = charFuncArgs[0].multiply(beta).add(kappa2).multiply(-2);
		Complex L34 = totalInstVar.subtract(kappa1 + kappa2 * theta)
				.subtract(charFuncArgs[0].multiply(2 * theta * beta)).multiply(3);
		Complex L35 = mixedDeg3.multiply(3).subtract(charFuncArgs[0].multiply(Math.pow(theta, 2) * beta)).multiply(4);
		Complex[] L3 = { L31, L32, L33, L34, L35 };
		// L[4, 3], L[4, 4] = -3.0 * (kappa2_p + beta * vol_backbone_eta * phi),
		// 2.0 * (vartheta2 - 2.0 * kappa_p - 4.0 * theta * beta * vol_backbone_eta * phi)*/
		Complex L41 = complex0;
		Complex L42 = complex0;
		Complex L43 = complex0;
		Complex L44 = charFuncArgs[0].multiply(beta).add(kappa2).multiply(-3);
		// TODO: In next line: Should be mulitplied with 3 instead of one
		Complex L45 = totalInstVar.multiply(1).subtract(2 * (kappa1 + kappa2 * theta))
				.subtract(charFuncArgs[0].multiply(4 * theta * beta)).multiply(2);
		Complex[] L4 = { L41, L42, L43, L44, L45 };

		// 3. Scalars
		Complex H0 = charFuncArgs[0].multiply(charFuncArgs[0]).add(charFuncArgs[0])
				.subtract(charFuncArgs[1].multiply(2)).multiply(0.5 * Math.pow(theta, 2));
		Complex H1 = charFuncArgs[0].multiply(charFuncArgs[0]).add(charFuncArgs[0])
				.subtract(charFuncArgs[1].multiply(2)).multiply(theta);
		Complex H2 = charFuncArgs[0].multiply(charFuncArgs[0]).add(charFuncArgs[0])
				.subtract(charFuncArgs[1].multiply(2)).multiply(0.5);
		Complex H3 = complex0;
		Complex H4 = complex0;

		// 2. Functions
		BiFunction<Double, Complex[], Complex> A0 = new BiFunction<Double, Complex[], Complex>() {

			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
                /*Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M0, complexes))
                        .add(LNSVQDUtils.scalarProduct(L0, complexes)).add(H0);*/
				Complex result = M0[1][1].multiply(complexes[1]).multiply(complexes[1])
						.add(L0[1].multiply(complexes[1])).add(L0[2].multiply(complexes[2])).add(H0);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A1 = new BiFunction<Double, Complex[], Complex>() {

			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
                /*Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M1, complexes))
                        .add(LNSVQDUtils.scalarProduct(L1, complexes)).add(H1);*/
				Complex result = M1[1][1].multiply(complexes[1]).multiply(complexes[1])
						.add(M1[1][2].multiply(complexes[1]).multiply(complexes[2]).multiply(2))
						.add(L1[1].multiply(complexes[1])).add(L1[2].multiply(complexes[2]))
						.add(L1[3].multiply(complexes[3])).add(H1);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A2 = new BiFunction<Double, Complex[], Complex>() {

			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
                /*Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M2, complexes))
                        .add(LNSVQDUtils.scalarProduct(L2, complexes)).add(H2);*/
				Complex result = M2[1][1].multiply(complexes[1]).multiply(complexes[1])
						.add(M2[2][2].multiply(complexes[2]).multiply(complexes[2]))
						.add(M2[1][2].multiply(complexes[1]).multiply(complexes[2]).multiply(2))
						.add(M2[1][3].multiply(complexes[1]).multiply(complexes[3]).multiply(2))
						.add(L2[1].multiply(complexes[1])).add(L2[2].multiply(complexes[2]))
						.add(L2[3].multiply(complexes[3])).add(L2[4].multiply(complexes[4])).add(H2);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A3 = new BiFunction<Double, Complex[], Complex>() {

			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
                /*Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M3, complexes))
                        .add(LNSVQDUtils.scalarProduct(L3, complexes)).add(H3);*/
				Complex result = M3[2][2].multiply(complexes[2]).multiply(complexes[2])
						.add(M3[1][2].multiply(complexes[1]).multiply(complexes[2]).multiply(2))
						.add(M3[1][3].multiply(complexes[1]).multiply(complexes[3]).multiply(2))
						.add(M3[1][4].multiply(complexes[1]).multiply(complexes[4]).multiply(2))
						.add(M3[2][3].multiply(complexes[2]).multiply(complexes[3]).multiply(2))
						.add(L3[2].multiply(complexes[2])).add(L3[3].multiply(complexes[3]))
						.add(L3[4].multiply(complexes[4])).add(H3);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A4 = new BiFunction<Double, Complex[], Complex>() {

			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
                /*Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M4, complexes))
                        .add(LNSVQDUtils.scalarProduct(L4, complexes)).add(H4);*/
				Complex result = M4[2][2].multiply(complexes[2]).multiply(complexes[2])
						.add(M4[3][3].multiply(complexes[3]).multiply(complexes[3]))
						.add(M4[1][3].multiply(complexes[1]).multiply(complexes[3]).multiply(2))
						.add(M4[1][4].multiply(complexes[1]).multiply(complexes[4]).multiply(2))
						.add(M4[2][3].multiply(complexes[2]).multiply(complexes[3]).multiply(2))
						.add(M4[2][4].multiply(complexes[2]).multiply(complexes[4]).multiply(2))
						.add(L4[3].multiply(complexes[3])).add(L4[4].multiply(complexes[4])).add(H4);
				return result;
			}
		};
		ArrayList<BiFunction<Double, Complex[], Complex>> As = new ArrayList<>();
		As.add(0, A0);
		As.add(1, A1);
		As.add(2, A2);
		As.add(3, A3);
		As.add(4, A4);
		return As;
	}

	// The next method return the As solution path
	public Complex[][] getSolutionPathForODESystem(double[] timeGridForPathIntegration, Complex[] charFuncArgs) {
		/**
		 * The ODE-system for the second-order exponential-affine approximation
		 */
		List<BiFunction<Double, Complex[], Complex>> odeSystem = getOdeSystem(charFuncArgs);

		/**
		 * Calculate the solution for the A's
		 */
		Complex[] state = new Complex[]{Complex.ZERO, charFuncArgs[2].multiply(-1), Complex.ZERO, Complex.ZERO, Complex.ZERO};
		// Solve ODE for A^{(k)}'s
		ComplexRungeKutta4thOrderIntegrator complexRungeKutta4thOrderIntegrator = new ComplexRungeKutta4thOrderIntegrator(state, odeSystem);
		Complex[][] solutionPath = complexRungeKutta4thOrderIntegrator.getSolutionPath(timeGridForPathIntegration);

		return solutionPath;
	}

	// Calculate the affine-exponential approximation to the characteristic function (whole path)
	public Complex[] calculateExponentialAffineApproximationFullPath(double[] timeGridForPathIntegration, Complex[] charFuncArgs) {
		// Fetch the As
		Complex[][] A = getSolutionPathForODESystem(timeGridForPathIntegration, charFuncArgs);

		// Get the approximation to the MGF for every time-point in our time-grid
		Complex[] complexAffineApproximationPath = new Complex[timeGridForPathIntegration.length];
		for(int j = 0; j < timeGridForPathIntegration.length; j++) {
			Complex result = (charFuncArgs[0].multiply(-X0)
					.add(charFuncArgs[1].multiply(-I0))
					.add(A[j][0])
					.add(A[j][1].multiply(Y0))
					.add(A[j][2].multiply(Math.pow(Y0, 2)))
					.add(A[j][3].multiply(Math.pow(Y0, 3)))
					.add(A[j][4].multiply(Math.pow(Y0, 4))))
					.exp();
			// TODO: Check and put somewhere else
			complexAffineApproximationPath[j] = result.isNaN() ? Complex.ZERO : result;
		}
		return complexAffineApproximationPath;
	}

	// Calculate the affine-exponential approximation to the characteristic function
	public Complex calculateExponentialAffineApproximation(Double endTime, Complex[] charFuncArgs) {
		int numStepsForODEIntegration = (int) (endTime * 365 * numStepsForODEIntegrationPerYear);
		// double
		double[] timeGridForMGFApproximationCalculation = LNSVQDUtils.createTimeGrid(0, endTime, numStepsForODEIntegration);
		// Choose the end point of the solution path
		Complex result = calculateExponentialAffineApproximationFullPath(timeGridForMGFApproximationCalculation, charFuncArgs)[numStepsForODEIntegration];
		return result;
	}

	/**
	 * We use our Runge-Kutta implementation to calculate the integral
	 */
	public double getCallPrice(double strike, double ttm) throws Exception {
		List<Pair<Double, Double>> strikeMaturityPairs = new ArrayList<>();
		strikeMaturityPairs.add(new Pair<>(ttm, strike));
		return getEuropeanOptionPrices(strikeMaturityPairs, true)[0];
	}

	public double[] getEuropeanOptionPrices(List<Pair<Double, Double>> strikeMaturityPairs, Boolean isCall) throws Exception {
		double[] uS = getU(strikeMaturityPairs);
		double[] optionPrices;

		if(isCall) {
			optionPrices = Arrays.stream(uS).map(u -> spot0 - u).toArray();
		} else {
			optionPrices = IntStream.range(0, strikeMaturityPairs.size())
					.mapToDouble(i -> {
						double maturity = strikeMaturityPairs.get(i).getKey();
						double strike = strikeMaturityPairs.get(i).getValue();
						double discountFactor = equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
						return discountFactor * strike - uS[i];
					})
					.toArray();
		}
		return optionPrices;
	}

	public double[] getU(List<Pair<Double, Double>> strikeMaturityPairs) throws Exception {
		int numStepsForODEIntegration = (int) (Math.max(strikeMaturityPairs.get(strikeMaturityPairs.size() - 1).getKey() * numStepsForODEIntegrationPerYear, 300));

		// Initialize the array of option prices that will be returned
		double[] optionPrices = new double[strikeMaturityPairs.size()];

		/**
		 * 1a. Extract time information from strike-maturity pairs
		 */
		double[] maturitiesWithZero = Stream.concat(strikeMaturityPairs.stream(), Stream.of(new Pair<Double, Double>(0., 0.)))
				.mapToDouble(pair -> pair.getFirst()).sorted().distinct().toArray();
		List<Double> maturitiesWithZeroList = Arrays.stream(maturitiesWithZero).boxed().collect(Collectors.toList());

		// Add points
		List<Double> timeGridForMGFApproximationCalculationList = Arrays.equals(maturitiesWithZero, new double[]{0}) ?
				Arrays.stream(maturitiesWithZero).boxed().collect(Collectors.toList()) :
				LNSVQDUtils.addTimePointsToArray(maturitiesWithZero, numStepsForODEIntegration + 1 - maturitiesWithZero.length);
		double[] timeGridForMGFApproximationCalculation = timeGridForMGFApproximationCalculationList.stream().mapToDouble(Double::doubleValue).toArray();

		/**
		 * 1b. Extract information for ODE-integration
		 */
		List<Double> gridForIntegrationWithMidPoints = LNSVQDUtils.addMidPointsToList(yGridForIntegration);

		/**
		 * 2. Precalcuate the expAffApprox path
		 * IMPORTANT: We need to calculate for all realizations of charFuncArgs!
		 *
		 * Result: The first dimension refers to a single realization of charFuncArgs, the second one to the value at the corresponding time
		 *
		 * NOTE: expAffApproxMatPathPerCharFuncRealization has the points yGridForInfiniteIntegral + the set of all midpoints because of how
		 * Runge-Kutta is implemented; hence, we have to evaluate E2 on the midpoints, too, i.e. on the points conatined in the array
		 * gridForIntegrationWithMidPoints
		 */
		// expAffApproxMatPathPerCharFuncArgRealization = Exponential-affine approximation at maturities per realization of the characteristic functions args
		Complex[][] expAffApproxMatPathPerCharFuncRealization = new Complex[gridForIntegrationWithMidPoints.size()][maturitiesWithZero.length];
		// Loop over charFunc args
		for(int l = 0; l < gridForIntegrationWithMidPoints.size(); l++) {
			// TODO: expAffApproxPathPerCharFuncRealization[l] should only contain values for maturity-points, not for filler points!
			double y = gridForIntegrationWithMidPoints.get(l);
			Complex[] charFuncArs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
			Complex[] expAffApproxPathPerCharFuncRealization = calculateExponentialAffineApproximationFullPath(timeGridForMGFApproximationCalculation, charFuncArs);

			// Get only maturities
			for(int i = 0; i < maturitiesWithZero.length; i++) {
				double time = maturitiesWithZero[i];
				int maturityIndex = timeGridForMGFApproximationCalculationList.indexOf(time);
				expAffApproxMatPathPerCharFuncRealization[l][i] = expAffApproxPathPerCharFuncRealization[maturityIndex];
			}
		}

		/**
		 * 3. For every maturity and every strike, we calculate the call option price
		 */
		for(int i = 0; i < strikeMaturityPairs.size(); i++) {
			double maturity = strikeMaturityPairs.get(i).getKey();
			double strike = strikeMaturityPairs.get(i).getValue();

			double forward = equityForwardStructure.getForward(maturity);
			double discountFactor = equityForwardStructure.getRepoCurve().getDiscountFactor(maturity); //Math.exp(-getRiskFreeRate(maturity) * maturity);
			double convenienceFactor = 0;

			// Get the time index of the maturity
			int maturityIndex = maturitiesWithZeroList.indexOf(maturity);
			int optionIndex = i;

			double logMoneyness = Math.log(spot0 / strike) + (-Math.log(discountFactor) - convenienceFactor);

			/**
			 * TODO: Change factor
			 * Define the function y -> factor * E2(t, -0.5 + iy)
			 */
			DoubleUnaryOperator integrand = new DoubleUnaryOperator() {
				@Override
				public double applyAsDouble(double y) {
					int yIndex = gridForIntegrationWithMidPoints.indexOf(y);
					if(yIndex == -1) {
						throw new IllegalArgumentException("y not in array: The E2 value for this the y-value " + y + " hasn't been calculated!");
					}

					// 1. Compute the value of the affine-exponential approximation
					Complex E2 = expAffApproxMatPathPerCharFuncRealization[yIndex][maturityIndex].multiply(new Complex(-0.5, y).multiply(X0).exp());

					// 2. Calculate result
					Complex result = new Complex(0.5, -y).multiply(logMoneyness).exp().multiply(1 / (y * y + 0.25)).multiply(E2);

					// 3. Get real part and return
					double resultReal = result.getReal();
					return resultReal;
				}
			};
			/**
			 * Integerate
			 */
			double result = discountFactor * (strike / Math.PI) * simpsonRealIntegrator.integrate(integrand);
			optionPrices[optionIndex] = result;

			/**
			 * *********************
			 * END
			 * **********************
			 */

		}
		/**
		 * 4. Return prices
		 */
		return optionPrices;
	}

	/**
	 * ***************************************************+
	 * SECTION 2: Get implied vol surface
	 * ***************************************************+
	 */
	public double[] getImpliedVolsStrikeMatList(ArrayList<Pair<Double, Double>> strikeMaturityPairs, double[] maturities) throws Exception {
		double[] impliedVols = new double[strikeMaturityPairs.size()];

		// Check if ...
		double[] uS = getU(strikeMaturityPairs);

		for(int i = 0; i < uS.length; i++) {
			double strike = strikeMaturityPairs.get(i).getValue();
			double ttm = strikeMaturityPairs.get(i).getKey();
			double discountFactor = equityForwardStructure.getRepoCurve().getDiscountFactor(ttm);
			double forward = equityForwardStructure.getForward(ttm);

			double impliedVol;
			double u = uS[i];
			// Use put-prices for itm-call region
			if(strike > forward) {
				double price = spot0 - u;
				impliedVol = Black76Model.optionImpliedVolatility(forward, strike, ttm, price / discountFactor, true);
			} else {
				double price = discountFactor * strike - u;
				impliedVol = Black76Model.optionImpliedVolatility(forward, strike, ttm, price / discountFactor, false);
			}
			impliedVols[i] = impliedVol;
		}
		return impliedVols;
	}

	// Method takes an existing volatility surface and creates a model implied vol-surface with the same structure
	public VolatilityPointsSurface getImpliedVolSurfaceFromVolSurface(VolatilityPointsSurface volatilityPointsSurface, double[] maturities) throws Exception {
		// Check if ...
		assert (volatilityPointsSurface.getToday() == spotDate) : "Expected spotDate " + spotDate + ", but got " + volatilityPointsSurface.getToday();;

		LocalDate today = volatilityPointsSurface.getToday();
		DayCountConvention dayCountConvention = volatilityPointsSurface.getDayCountConvention();

		// Get strike maturity pairs
		ArrayList<Pair<Double, Double>> strikeMaturityPairs = new ArrayList<>();
		for(VolatilityPoint volatilityPoint : volatilityPointsSurface.getVolatilityPoints()) {
			LocalDate maturity = volatilityPoint.getDate();
			double strike = volatilityPoint.getStrike();
			double ttm = dayCountConvention.getDaycountFraction(today, maturity);
			strikeMaturityPairs.add(new Pair<>(ttm, strike));
		}

		double[] uS = getU(strikeMaturityPairs);

		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();
		for(int i = 0; i < uS.length; i++) {
			LocalDate maturity = volatilityPointsSurface.getVolatilityPoints().get(i).getDate();
			double strike = volatilityPointsSurface.getVolatilityPoints().get(i).getStrike();
			double ttm = dayCountConvention.getDaycountFraction(today, maturity);
			double discountFactor = equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
			double forward = equityForwardStructure.getForward(ttm);

			double impliedVol;
			double u = uS[i];
			// Use put-prices for itm-call region
			if(strike > forward) {
				double price = spot0 - u;
				impliedVol = Black76Model.optionImpliedVolatility(forward, strike, ttm, price / discountFactor, true);
			} else {
				double price = discountFactor * strike - u;
				impliedVol = Black76Model.optionImpliedVolatility(forward, strike, ttm, price / discountFactor, false);
			}
			volatilityPoints.add(new VolatilityPoint(maturity, strike, impliedVol));
		}

		return new VolatilityPointsSurface(volatilityPoints, today, dayCountConvention);
	}

	public void setYGridForIntegration(double lb, double ub, int numberOfEvaluationPoints) {
		/**
		 * Get all the integration points from the integrator, in our case Simpson
		 */
		// Next lines adapted from finmath's Simpson implementation
		final double range = ub - lb;

		final int numberOfDoubleSizeIntervals = (int) ((numberOfEvaluationPoints - 1) / 2.0);

		final double doubleInterval = range / numberOfDoubleSizeIntervals;
		final double singleInterval = 0.5 * doubleInterval;

		IntStream intervals = IntStream.range(1, numberOfDoubleSizeIntervals);

		intervals.forEach(
				i -> {
					yGridForIntegration.add(lb + i * doubleInterval);
					yGridForIntegration.add(lb + i * doubleInterval + singleInterval);
				}
		);

		yGridForIntegration.add(lb + singleInterval);
		yGridForIntegration.add(lb);
		yGridForIntegration.add(ub);

		yGridForIntegration = yGridForIntegration.stream().sorted().distinct().collect(Collectors.toList());

	}

}
