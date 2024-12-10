package net.finmath.equities.models;

import org.apache.commons.math3.complex.Complex;

import java.util.HashMap;
import java.util.List;
import java.util.function.BiFunction;

public class ComplexRungeKutta4thOrderIntegrator {

	private Complex[] state;
	private List<BiFunction<Double, Complex[], Complex>> odeSystem;

	public ComplexRungeKutta4thOrderIntegrator(Complex[] state, List<BiFunction<Double, Complex[], Complex>> odeSystem) {
		this.state = state;
		this.odeSystem = odeSystem;
	}

	public void setState(Complex[] state) {
		this.state = state;
	}

	public Complex[] getState() {
		return state;
	}

	public Complex[][] getSolutionPath(double[] timePoints) {
		Complex[][] solutionPath = new Complex[timePoints.length][state.length];
		solutionPath[0] = this.state;
		for(int j = 1; j < timePoints.length; j++) {
			solutionPath[j] = integrateOneStep(timePoints[j - 1], timePoints[j] - timePoints[j - 1]);
			setState(solutionPath[j]);
		}
		return solutionPath;
	}

	private Complex[] integrateOneStep(double time, double timeStep) {
		Complex[] newState = this.state.clone();
		for(int j = 0; j < state.length; j++) {
			Complex[] overwrittenState = newState.clone();
			Complex d1 = this.odeSystem.get(j).apply(time, overwrittenState);

			overwrittenState[j] = newState[j].add(d1.multiply(timeStep / 2));
			Complex d2 = this.odeSystem.get(j).apply(time + (timeStep / 2), overwrittenState);

			overwrittenState[j] = newState[j].add(d2.multiply(timeStep / 2));
			Complex d3 = this.odeSystem.get(j).apply(time + (timeStep / 2), overwrittenState);

			overwrittenState[j] = newState[j].add(d3.multiply(timeStep / 2));
			Complex d4 = this.odeSystem.get(j).apply(time + timeStep, overwrittenState);

			newState[j] = newState[j].add((d1.add(d2.multiply(2)).add(d3.multiply(2)).add(d4).multiply(timeStep / 6)));
		}
		return newState;
	}

}
