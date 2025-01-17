package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.ComplexRungeKutta4thOrderIntegrator;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.RungeKuttaIntegrator;
import org.junit.Assert;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;

class ComplexRungeKutta4thOrderIntegratorTest {


	@Test
	void setState() {
	}

	@Test
	void getSolutionPath() {
		// Create ComplexRungeKutta4thOrderIntegrator for complex harmonic oscillation equation
		BiFunction<Double, Complex[], Complex> ode = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = complexes[0].multiply(-0.5).add(complexes[0].multiply(new Complex(0., 1.)).multiply(2));
				return result;
			}
		};

		final Complex[] state = {new Complex(1.0, 0.)};
		List<BiFunction<Double, Complex[], Complex>> odeSystem = new ArrayList<>();
		odeSystem.add(ode);
		ComplexRungeKutta4thOrderIntegrator odeRungeKutta = new ComplexRungeKutta4thOrderIntegrator(state, odeSystem);

		// TEST; TO DELETE
		//odeRungeKutta.setState(new Complex[]{new Complex(5.0, 0.)});
		//System.out.println("State: " + state[0].getReal());

		double t0 = 0.0;
		double stepSize = 0.01;
		int steps = 999;
		double[] timePoints = createTimeGrid(t0, stepSize, steps);

		Complex[][] solutionPath = odeRungeKutta.getSolutionPath(timePoints);

		//Analytical solution
		Function<Double, Complex> harmonicOscillation = new Function<Double, Complex>() {
			@Override
			public Complex apply(Double aDouble) {
				Complex r = new Complex(-0.5, 2);
				return state[0].multiply(r.multiply(aDouble).exp());
			}
		};

		// Print results
		System.out.println("Time" + "\t"+ "RungeKutta" + "\t\t" + "Analytical");
		System.out.println("Time" + "\t"+ "Real" + "\t" + "Imagniary"+ "\t"+ "Real" + "\t" + "Imagniary");
		for(int k = 0; k <= steps; k++) {
			double currentTime = timePoints[k];
			System.out.println(currentTime + "\t"+ solutionPath[k][0].getReal() + "\t" + solutionPath[k][0].getImaginary() + "\t" + harmonicOscillation.apply(currentTime).getReal() + "\t" + harmonicOscillation.apply(currentTime).getImaginary());
		}

		//TODO: Assert that relative error is below a certain threshold
	}

	@Test
	void testLotkaVolterra() {
		double t0 = 0.0;
		double tEnd = 0.2;
		int steps = 20;
		double stepSize = (tEnd - t0) / steps;
		double[] timePoints = LNSVQDUtils.createTimeGrid(t0, tEnd, steps);


		double alpha = 1;
		double beta = 0.1;
		double delta = 0.075;
		double gamma = 1.5;

		/**
		 * ***************************************************+
		 * 1. RK own implementation
		 * ***************************************************+
		 */
		BiFunction<Double, Complex[], Complex> evolutionPrey = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = complexes[0].multiply(alpha).subtract(complexes[0].multiply(complexes[1]).multiply(beta));
				return result;
			}
		};
		BiFunction<Double, Complex[], Complex> evolutionPredator = new BiFunction<Double, Complex[], Complex>() {
			@Override
			public Complex apply(Double aDouble, Complex[] complexes) {
				Complex result = complexes[0].multiply(complexes[1]).multiply(delta).subtract(complexes[1].multiply(gamma));
				return result;
			}
		};

		final Complex[] state = {new Complex(100, 0.), new Complex(100, 0.)};
		List<BiFunction<Double, Complex[], Complex>> odeSystem = new ArrayList<>();
		odeSystem.add(evolutionPrey);
		odeSystem.add(evolutionPredator);
		ComplexRungeKutta4thOrderIntegrator odeRungeKutta = new ComplexRungeKutta4thOrderIntegrator(state, odeSystem);

		Complex[][] solutionPath = odeRungeKutta.getSolutionPath(timePoints);

		// Print results
		System.out.println("Time" + "\t"+ "Prey" + "\t\t" + "Predator");
		System.out.println("Time" + "\t"+ "Real" + "\t"+ "Real");
		for(int k = 0; k <= steps; k++) {
			double currentTime = timePoints[k];
			System.out.println(currentTime + "\t"+ solutionPath[k][0].getReal() + "\t" + solutionPath[k][1].getReal());
		}

		/**
		 * ***************************************************+
		 * 2. RK Apache
		 * ***************************************************+
		 */
		FirstOrderDifferentialEquations firstOrderDifferentialEquations = new FirstOrderDifferentialEquations() {
			@Override
			public int getDimension() {
				return 2;
			}

			@Override
			public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
				yDot[0] = y[0] * (alpha) - (y[0] * y[1] * (beta));
				yDot[1] = (y[0] * (y[1]) * (delta)) - (y[1] * (gamma));
			}
		};

		double[] initalStateReal = LNSVQDUtils.convertComplexArrayToDoubleWithReal(state);
		double[] yEnd = initalStateReal.clone();
		ClassicalRungeKuttaIntegrator classicalRungeKuttaIntegrator = new ClassicalRungeKuttaIntegrator(stepSize);
		classicalRungeKuttaIntegrator.integrate(firstOrderDifferentialEquations, t0, initalStateReal, tEnd, yEnd);

		// Print results
		System.out.println("Final state at t = " + tEnd);
		System.out.println("y[0] (position): " + yEnd[0]);
		System.out.println("y[1] (velocity): " + yEnd[1]);


		/**
		 * ***************************************************+
		 * 3. Compare
		 * ***************************************************+
		 */
		Assertions.assertEquals(yEnd[0], solutionPath[timePoints.length - 1][0].getReal(), 1e-6);
		Assertions.assertEquals(yEnd[1], solutionPath[timePoints.length - 1][1].getReal(), 1e-6);
	}

	private double[] createTimeGrid(double t0, double stepSize, int steps) {
		double[] timeGrid = new double[steps + 1];
		timeGrid[0] = t0;
		for(int j = 1; j <= steps; j++) {
			timeGrid[j] = timeGrid[j - 1] + stepSize;
		}
		return timeGrid;
	}
}