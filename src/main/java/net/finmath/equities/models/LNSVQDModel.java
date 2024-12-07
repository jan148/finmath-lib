package net.finmath.equities.models;

import net.finmath.fouriermethod.CharacteristicFunction;
import net.finmath.integration.RealIntegral;
import net.finmath.marketdata.model.curves.DiscountCurve;
import org.apache.commons.math3.Field;
import org.apache.commons.math3.FieldElement;
import org.apache.commons.math3.complex.Complex;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaFieldIntegrator;
import org.apache.commons.math3.ode.nonstiff.RungeKuttaFieldIntegrator;
import org.apache.commons.math3.complex.ComplexField;

public class LNSVQDModel {
	/**
	 * Numerical parameters
	 */
	private final int numStepsForODEIntegration = 1000;
	private final long numStepsForInfiniteIntegral = 1000000000;
	private final double upperBoundForInfiniteIntegral = numStepsForInfiniteIntegral / 10;

	/**
	 * Model parameters under the EMM
	 */
	private final double spot0;
	private final double sigma0;
	private final double kappa1;
	private final double kappa2;
	private final double theta;
	private final double beta;
	private final double epsilon;

	/**
	 * Transformed inital values
	 */
	double X0, Y0, I0;

	public LNSVQDModel(double spot0, double sigma0, double kappa1, double kappa2, double theta, double beta, double epsilon, double I0) {
		this.spot0 = spot0;
		this.sigma0 = sigma0;
		this.kappa1 = kappa1;
		this.kappa2 = kappa2;
		this.theta = theta;
		this.beta = beta;
		this.epsilon = epsilon;

		this.X0 = Math.log(this.spot0);
		this.Y0 = this.sigma0 - this.theta;
		this.I0 = I0;

	}

	// Calculate the affine-exponential approximation to the characteristic function
	private Complex calculateExponentialAffineApproximation(Double ttm, Complex[] charFuncArgs){
		/**
		 * Create the functions for the exponential-affine approximation; Names like in the paper
		 */
		// 1. Matrices
		Complex[][] M0 = {{}, {}, {}, {}, {}};
		Complex[][] M1 = {{}, {}, {}, {}, {}};
		Complex[][] M2 = {{}, {}, {}, {}, {}};
		Complex[][] M3 = {{}, {}, {}, {}, {}};
		Complex[][] M4 = {{}, {}, {}, {}, {}};

		Complex[] L0 = {};
		Complex[] L1 = {};
		Complex[] L2 = {};
		Complex[] L3 = {};
		Complex[] L4 = {};

		Complex H0 = new Complex(0., 0.);
		Complex H1 = new Complex(0., 0.);
		Complex H2 = new Complex(0., 0.);
		Complex H3 = new Complex(0., 0.);
		Complex H4 = new Complex(0., 0.);

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
		Complex[] A = complexRungeKutta4thOrderIntegrator.getSolutionPath(LNSVQDUtils.createTimeGrid(0., ttm, this.numStepsForODEIntegration))[this.numStepsForODEIntegration];

		/**
		 * Calculate the solution for the A's
		 */
		Complex result = (charFuncArgs[0].multiply(-this.X0)
				.add(charFuncArgs[1].multiply(-this.I0))
				.add(A[0])
				.add(A[1].multiply(this.Y0))
				.add(A[2].multiply(Math.pow(this.Y0, 2)))
				.add(A[3].multiply(Math.pow(this.Y0, 3)))
				.add(A[4].multiply(Math.pow(this.Y0, 4))))
				.exp();

		return result;
	}

	/**
	 * We use our Runge-Kutta implementation to calculate the integral
	 */
	public double getCallPrice(double strike, double ttm, EquityForwardStructure equityForwardStructure, RealIntegral integrator){
		double discountFactor = equityForwardStructure.getGrowthDiscountFactor(0, ttm);
		double convenienceFactor = 0; //TODO: Replace by equityForwardStructure...
		double logMoneyness = Math.log(spot0 / strike) + Math.log(discountFactor - convenienceFactor);

		/**
		 * We use our Runge-Kutta implementation to calculate the integral
		 */
		BiFunction<Double, Complex[], Complex> integrand = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				// 1. Compute the value of the affine-exponential approximation
				Complex approxCharFunVal = calculateExponentialAffineApproximation(ttm, complexes);
				Complex E2 = approxCharFunVal
						.multiply(complexes[0].multiply(X0).add(complexes[1].multiply(I0)));

				// 2. Calcuate result
				Complex result = (new Complex(-0.5 * logMoneyness, aDouble * logMoneyness)).exp().multiply(1 / (aDouble * aDouble + 0.25)).multiply(E2);
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
		Complex integral = complexRungeKutta4thOrderIntegrator.getSolutionPath(LNSVQDUtils.createTimeGrid(0., upperBoundForInfiniteIntegral, this.numStepsForInfiniteIntegral))[this.numStepsForInfiniteIntegral][0];

		/**
		 * Get real part, multiply with factor
		 */
		double integralReal = integral.getReal();
		double optionPrice =  this.spot0 - (discountFactor * strike / Math.PI)  * integralReal;

		return optionPrice;
	}

}
