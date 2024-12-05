package net.finmath.equities.models;

import net.finmath.fouriermethod.CharacteristicFunction;
import net.finmath.integration.RealIntegral;
import net.finmath.marketdata.model.curves.DiscountCurve;
import org.apache.commons.math3.Field;
import org.apache.commons.math3.FieldElement;
import org.apache.commons.math3.complex.Complex;

import java.util.List;
import java.util.function.BiFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaFieldIntegrator;
import org.apache.commons.math3.ode.nonstiff.RungeKuttaFieldIntegrator;
import org.apache.commons.math3.complex.ComplexField;

public class LNSVQDModel {

	/**
	 * Model parameters under the EMM
	 */
	private final double spot0;
	private final double sigma0;
	private final double kappa1;
	private final double kappa2;
	private final double tehta;
	private final double beta;
	private final double epsilon;

	/**
	 * The ODE-system for the second-order exponential-affine approximation
	 */
	private final List<BiFunction<Double, Complex[], Complex>> odeSystem;

	/*

	 */
	private Function<Double[], Complex> exponentialAffineApproximation;

	public LNSVQDModel(double spot0, double sigma0, double kappa1, double kappa2, double tehta, double beta, double epsilon) {
		this.spot0 = spot0;
		this.sigma0 = sigma0;
		this.kappa1 = kappa1;
		this.kappa2 = kappa2;
		this.tehta = tehta;
		this.beta = beta;
		this.epsilon = epsilon;

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

		Complex L0 = ;
		Complex L1 = ;
		Complex L2 = ;
		Complex L3 = ;
		Complex L4 = ;

		// 2. Functions
		BiFunction<Double, Complex[], Complex> A0 = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex results = new Complex();
			}
		}

	}

	// Calculate the affine-exponential approximation to the characteristic function
	/*private void calculateExponentialAffineApproximation(){
		// 1. Describe the field

		// Solve PDE
		ClassicalRungeKuttaFieldIntegrator classicalRungeKuttaFieldIntegrator = new ClassicalRungeKuttaFieldIntegrator(ComplexField.getInstance(), new eTRealField(1.));
		// Function<Double, Comp>

		this.exponentialAffineApproximation = new Function<Double[], Complex>() {
			@Override
			public Complex apply(Double[] doubles) {
				return null;
			}
		};
	}*/

	/*private double getOptionPrice(double strike, double maturity, DiscountCurve discountCurve, EquityForwardStructure convenienceCurve, RealIntegral integrator){


		double discountFactor = discountCurve.getDiscountFactor(maturity);
		double convenienceFactor = convenienceCurve.getConvenienceFactor(maturity);
		double logMoneyness = Math.log(spot0 / strike) + Math.log(discountFactor - convenienceFactor);

		// Define the integrand
		DoubleUnaryOperator integrand = new DoubleUnaryOperator() {
			@Override
			public double applyAsDouble(double operand) {
				//TODO: Check next line; Among other things, the exp-affine approximation has to be replaced by E^{(2)}
				return (new Complex(-0.5 * logMoneyness, operand * logMoneyness)).exp().multiply(1 / (operand * operand + 0.25)).multiply(exponentialAffineApproximation.apply(new Double[]{0., 0., operand})).getReal();
			}

			@Override
			public DoubleUnaryOperator compose(DoubleUnaryOperator before) {
				return DoubleUnaryOperator.super.compose(before);
			}

			@Override
			public DoubleUnaryOperator andThen(DoubleUnaryOperator after) {
				return DoubleUnaryOperator.super.andThen(after);
			}
		};

		double integral = integrator.integrate(integrand);
		double optionPrice =  (discountFactor * strike / Math.PI)  * integral;

		return optionPrice;
	}*/

}
