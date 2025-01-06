package net.finmath.equities.models.LNSVQD;

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

		int n = newState.length;
		Complex[] d1 = new Complex[n];
		Complex[] d2 = new Complex[n];
		Complex[] d3 = new Complex[n];
		Complex[] d4 = new Complex[n];

		Complex[] overwrittenState = this.state.clone();
		for(int j = 0; j < n; j++) {
			d1[j] = odeSystem.get(j).apply(time, overwrittenState);
		}

		overwrittenState = this.state.clone();
		for(int j = 0; j < n; j++) {
			overwrittenState[j] = overwrittenState[j].add(d1[j].multiply(timeStep / 2));
		}
		for(int j = 0; j < n; j++) {
			d2[j] = odeSystem.get(j).apply(time + (timeStep / 2), overwrittenState);
		}

		overwrittenState = this.state.clone();
		for(int j = 0; j < n; j++) {
			overwrittenState[j] = overwrittenState[j].add(d2[j].multiply(timeStep / 2));
		}
		for(int j = 0; j < state.length; j++) {
			d3[j] = odeSystem.get(j).apply(time + (timeStep / 2), overwrittenState);
		}

		overwrittenState = this.state.clone();
		for(int j = 0; j < n; j++) {
			overwrittenState[j] = overwrittenState[j].add(d3[j].multiply(timeStep));
		}
		for(int j = 0; j < n; j++) {
			d4[j] = odeSystem.get(j).apply(time + timeStep, overwrittenState);
		}

		overwrittenState = this.state.clone();
		for(int j = 0; j < state.length; j++) {
			newState[j] = overwrittenState[j].add((d1[j].add(d2[j].multiply(2)).add(d3[j].multiply(2)).add(d4[j]).multiply(timeStep / 6)));
		}

		return newState;
	}

}
