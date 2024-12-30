package net.finmath.equities.models.LNSVQD;

import org.apache.commons.math3.complex.Complex;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

public class LNSVQDModelAnalyticalPricer extends LNSVQDModel {
	/**
	 * Numerical parameters
	 */
	private final int numStepsForODEIntegration = 100;
	private final int numStepsForInfiniteIntegral = 1000;
	private final double upperBoundForInfiniteIntegral = numStepsForInfiniteIntegral / 10;

	public LNSVQDModelAnalyticalPricer(double spot0, double sigma0, double kappa1, double kappa2, double theta, double beta, double epsilon, double I0) {
		super(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);
	}

	/**
	 * ***************************************************+
	 * SECTION 1: Semi-analytical call option price calculation
	 * ***************************************************+
	 */

	// Calculate the affine-exponential approximation to the characteristic function
	public Complex calculateExponentialAffineApproximation(Double ttm, Complex[] charFuncArgs){
		/**
		 * Create the functions for the exponential-affine approximation; Names like in the paper
		 */
		// Some constants
		Complex complex0 = Complex.ZERO;
		Complex totalInstVar = new Complex(this.totalInstVar, 0);
		Complex mixedDeg3 = new Complex(this.theta * this.totalInstVar, 0);
		Complex mixedDeg4 = new Complex(Math.pow(this.theta, 2) * this.totalInstVar, 0);

		// 1. Matrices
		Complex[][] M0 = {{complex0, complex0, complex0, complex0, complex0},
				{complex0, mixedDeg4.multiply(0.5), complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0}};

		Complex[][] M1 = {{complex0, complex0, complex0, complex0, complex0},
				{complex0, mixedDeg3, mixedDeg4, complex0, complex0},
				{complex0, mixedDeg4, complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0}};

		Complex[][] M2 = {{complex0, complex0, complex0, complex0, complex0},
				{complex0, totalInstVar.multiply(1. / 2), mixedDeg3.multiply(2), mixedDeg4.multiply(3. / 2), complex0},
				{complex0, mixedDeg3.multiply(2), mixedDeg4.multiply(2), complex0, complex0},
				{complex0, mixedDeg4.multiply(3. / 2), complex0, complex0, complex0},
				{complex0, complex0, complex0, complex0, complex0}};

		Complex[][] M3 = {{complex0, complex0, complex0, complex0, complex0},
				{complex0, complex0, totalInstVar, mixedDeg3.multiply(3), mixedDeg4.multiply(2)},
				{complex0, totalInstVar, mixedDeg3.multiply(4), mixedDeg4.multiply(3), complex0},
				{complex0, mixedDeg3.multiply(3), mixedDeg4.multiply(3), complex0, complex0},
				{complex0, mixedDeg4.multiply(2), complex0, complex0, complex0}};

		Complex[][] M4 = {{complex0, complex0, complex0, complex0, complex0},
				{complex0, complex0, complex0, totalInstVar.multiply(3. / 2), mixedDeg3.multiply(4)},
				{complex0, complex0, totalInstVar.multiply(2), mixedDeg3.multiply(6), mixedDeg4.multiply(4)},
				{complex0, totalInstVar.multiply(3. / 2), mixedDeg3.multiply(6), mixedDeg4.multiply(9. / 2), complex0},
				{complex0, mixedDeg3.multiply(4), mixedDeg4.multiply(4), complex0, complex0}};

		// 2. Vectors
		Complex L01 = complex0;
		Complex L02 = charFuncArgs[0].multiply(-Math.pow(theta, 2) * beta);
		Complex L03 = mixedDeg4;
		Complex L04 = complex0;
		Complex L05 = complex0;
		Complex[] L0 = {L01, L02, L03, L04, L05};

		Complex L11 = complex0;
		Complex L12 = charFuncArgs[0].multiply(-theta * beta * 2).subtract(kappa1 + kappa2 * theta);
		Complex L13 = mixedDeg3.subtract(charFuncArgs[0].multiply(Math.pow(theta, 2) * beta)).multiply(2);
		Complex L14 = mixedDeg4.multiply(3);
		Complex L15 = complex0;
		Complex[] L1 = {L11, L12, L13, L14, L15};

		Complex L21 = complex0;
		Complex L22 = charFuncArgs[0].multiply(-beta).subtract(kappa2);
		Complex L23 = totalInstVar.subtract(2 * (kappa1 + kappa2 * theta)).subtract(charFuncArgs[0].multiply(theta * beta).multiply(4));
		Complex L24 = totalInstVar.multiply(-2).subtract(charFuncArgs[0].multiply(Math.pow(theta, 2) * beta)).multiply(3); // q = -p = -1
		Complex L25 = mixedDeg4.multiply(6);
		Complex[] L2 = {L21, L22, L23, L24, L25};

		Complex L31 = complex0;
		Complex L32 = complex0;
		Complex L33 = charFuncArgs[0].multiply(beta).add(kappa2).multiply(-2);
		Complex L34 = totalInstVar.subtract(kappa1 + kappa2 * theta).subtract(charFuncArgs[0].multiply(2 * theta * beta)).multiply(3);
		Complex L35 = mixedDeg3.multiply(3).subtract(charFuncArgs[0].multiply(Math.pow(theta, 2) * beta)).multiply(4);
		Complex[] L3 = {L31, L32, L33, L34, L35};

		Complex L41 = complex0;
		Complex L42 = complex0;
		Complex L43 = complex0;
		Complex L44 = charFuncArgs[0].multiply(beta).add(kappa2).multiply(-3);
		Complex L45 = totalInstVar.multiply(3).subtract(2 * (kappa1 + kappa2 * theta)).subtract(charFuncArgs[0].multiply(4 * theta * beta)).multiply(2);
		Complex[] L4 = {L41, L42, L43, L44, L45};

		// 3. Scalars
		Complex H0 = charFuncArgs[0].pow(2).add(charFuncArgs[0]).subtract(charFuncArgs[1].multiply(2)).multiply(Math.pow(theta, 2) / 2);
		Complex H1 = charFuncArgs[0].pow(2).add(charFuncArgs[0]).subtract(charFuncArgs[1].multiply(2)).multiply(theta);
		Complex H2 = charFuncArgs[0].pow(2).add(charFuncArgs[0]).subtract(charFuncArgs[1].multiply(2)).multiply(0.5);
		Complex H3 = complex0;
		Complex H4 = complex0;

		// 2. Functions
		BiFunction<Double, Complex[], Complex> A0 = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M0, complexes)).add(LNSVQDUtils.scalarProduct(L0, complexes)).add(H0);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A1 = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M1, complexes)).add(LNSVQDUtils.scalarProduct(L1, complexes)).add(H1);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A2 = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M2, complexes)).add(LNSVQDUtils.scalarProduct(L2, complexes)).add(H2);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A3 = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M3, complexes)).add(LNSVQDUtils.scalarProduct(L3, complexes)).add(H3);
				return result;
			}
		};

		BiFunction<Double, Complex[], Complex> A4 = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = LNSVQDUtils.scalarProduct(complexes, LNSVQDUtils.matrixVectorMult(M4, complexes)).add(LNSVQDUtils.scalarProduct(L4, complexes)).add(H4);
				return result;
			}
		};

		/**
		 * The ODE-system for the second-order exponential-affine approximation
		 */
		List<BiFunction<Double, Complex[], Complex>> odeSystem = new ArrayList<>();
		odeSystem.add(A0);
		odeSystem.add(A1);
		odeSystem.add(A2);
		odeSystem.add(A3);
		odeSystem.add(A4);

		/**
		 * Calculate the solution for the A's
		 */
		Complex[] state = new Complex[]{new Complex(0., 0.), charFuncArgs[2].multiply(-1), new Complex(0., 0.), new Complex(0., 0.), new Complex(0., 0.)};
		// Solve ODE for A^{(k)}'s
		ComplexRungeKutta4thOrderIntegrator complexRungeKutta4thOrderIntegrator = new ComplexRungeKutta4thOrderIntegrator(state, odeSystem);
		// Choose the end point of the solution path
		Complex[] A = complexRungeKutta4thOrderIntegrator.getSolutionPath(LNSVQDUtils.createTimeGrid(0, ttm, this.numStepsForODEIntegration))[this.numStepsForODEIntegration];

		/**
		 * Calculate the solution for the A's
		 */
		Complex result = (charFuncArgs[0].multiply(-X0)
				.add(charFuncArgs[1].multiply(-I0))
				.add(A[0])
				.add(A[1].multiply(Y0))
				.add(A[2].multiply(Math.pow(Y0, 2)))
				.add(A[3].multiply(Math.pow(Y0, 3)))
				.add(A[4].multiply(Math.pow(Y0, 4))))
				.exp();

		return result;
	}

	/**
	 * We use our Runge-Kutta implementation to calculate the integral
	 */
	public double getCallPrice(double strike, double ttm, double discountFactor, double convencienceFactor){
		double logMoneyness = Math.log(spot0 / strike) + (-Math.log(discountFactor) - convencienceFactor);

		/**
		 * We use our Runge-Kutta implementation to calculate the integral
		 */
		BiFunction<Double, Complex[], Complex> integrand = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, aDouble), Complex.ZERO, Complex.ZERO};

				// 1. Compute the value of the affine-exponential approximation
				Complex approxCharFuncVal = calculateExponentialAffineApproximation(ttm, charFuncArgs);
				Complex E2 = approxCharFuncVal
						.multiply(charFuncArgs[0].multiply(X0).add(charFuncArgs[1].multiply(I0).exp()));

				// 2. Calculate result
				Complex result = new Complex(0.5, -aDouble).multiply(logMoneyness).exp().multiply(1 / (aDouble * aDouble + 0.25)).multiply(E2);
				return result;
			}
		};

		/**
		 * Integerate
		 */
		Complex[] state = new Complex[]{new Complex(0., 0.)};
		List<BiFunction<Double, Complex[], Complex>> odeSystem = new ArrayList<>();
		odeSystem.add(integrand);
		ComplexRungeKutta4thOrderIntegrator complexRungeKutta4thOrderIntegrator = new ComplexRungeKutta4thOrderIntegrator(state, odeSystem);
		// Choose the end point of the solution path
		Complex integral = complexRungeKutta4thOrderIntegrator.getSolutionPath(LNSVQDUtils.createTimeGrid(0, upperBoundForInfiniteIntegral, this.numStepsForInfiniteIntegral))[this.numStepsForInfiniteIntegral][0];

		/**
		 * Get real part, multiply with factor
		 */
		double integralReal = integral.getReal();
		double optionPrice =  spot0 - (discountFactor * strike / Math.PI) * integralReal;

		return optionPrice;
	}

}
