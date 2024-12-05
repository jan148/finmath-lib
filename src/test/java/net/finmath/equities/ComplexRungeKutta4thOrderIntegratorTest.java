package net.finmath.equities;

import net.finmath.equities.models.ComplexRungeKutta4thOrderIntegrator;
import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

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

	private double[] createTimeGrid(double t0, double stepSize, int steps) {
		double[] timeGrid = new double[steps + 1];
		timeGrid[0] = t0;
		for(int j = 1; j <= steps; j++) {
			timeGrid[j] = timeGrid[j - 1] + stepSize;
		}
		return timeGrid;
	}
}