package net.finmath.equities.models;

import org.apache.commons.math3.complex.Complex;

import java.util.List;
import java.util.function.BiFunction;

public class ComplexRungeKutta4thOrderIntegrator {

	Complex[] state;
	List<BiFunction<Double, Complex[], Complex>> odeSystem;

	public ComplexRungeKutta4thOrderIntegrator(List<BiFunction<Double, Complex[], Complex>> odeSystem) {
		this.state = state;
		this.odeSystem = odeSystem;
	}

	public void setState(Complex[] state) {
		this.state = state;
	}

	private Complex[] integrateOneStep(double time, double timeStep){
		Complex[] newState = new Complex[state.length];
		for(int j = 0; j < state.length; j++) {
			Complex[] overwrittenState = this.state;
			Complex d1 = this.odeSystem.get(j).apply(time, state);

			overwrittenState[j] = this.state[j].add(d1.multiply(timeStep / 2));
			Complex d2 = this.odeSystem.get(j).apply(time / 2, overwrittenState);

			overwrittenState[j] = this.state[j].add(d2.multiply(timeStep / 2));
			Complex d3 = this.odeSystem.get(j).apply(time / 2, overwrittenState);

			overwrittenState[j] = this.state[j].add(d3.multiply(timeStep / 2));
			Complex d4 = this.odeSystem.get(j).apply(time, overwrittenState);

			newState[j] = state[j].add((d1.add(d2.multiply(2)).add(d3.multiply(2)).add(d4)).multiply(timeStep / 6));
		}
		return newState;
	}



}
