package net.finmath.montecarlo;

import net.finmath.equities.models.LNSVQDModel;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.LNSVQDDiscretizationScheme;
import net.finmath.rootfinder.NewtonsMethod;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import org.junit.Assert;
import org.junit.Test;

import java.util.function.Function;

public class LNSVQDDiscretizationSchemeTest {
	/**
	 * Simulation parameters
	 */
	int numberOfPaths = 26;
	int seed = 3107;

	/**
	 * Model params
	 */
	private final double spot0 = 1;
	private final double sigma0 = 0.5;
	private final double kappa1 = 0.1;
	private final double kappa2 = 0.0;
	private final double theta = 0.5;
	private final double beta = 0.2;
	private final double epsilon = 0;


	/**
	 * Models
	 */
	LNSVQDModel lnsvqdModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);

	/**
	 * Option params
	 */
	double strike = 1;
	double maturity = 1. / 12;

	/**
	 * Market observables
	 */
	double riskFreeRate = 0.05;
	double discountFactor = Math.exp(-riskFreeRate * maturity);
	double convenienceFcator = 0;

	/**
	 * Time discretization
	 */
	double[] timeGrid = LNSVQDUtils.createTimeGrid(0, 5, 100);
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
	public void doPrecalculateProcess() throws CalculationException {
		// Create 2D-Brownian motion and discretization scheme
		BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
		LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdModel, brownianMotion);

		for(int l = 0; l < lnsvqdModel.getNumberOfComponents(); l++) {
			System.out.print("Time\t");
			for(int j = 0; j < numberOfPaths; j++) {
				System.out.print("Realization " + j + "\t");
			}
			for(int k = 0; k < timeDiscretization.getNumberOfTimes(); k++) {
				System.out.print("\n" + timeDiscretization.getTime(k) + "\t");
				for(int j = 0; j < numberOfPaths; j++) {
					System.out.print(lnsvqdDiscretizationScheme.getProcessValue(k, l).get(j) + "\t");
				}
			}
			System.out.print("\n\n");
		}
	}

	@Test
	public void printRootFunction() {
		double c = 10;
		double deltaT = 0.5;
		Function<RandomVariable, Double> rootFunction = new Function<RandomVariable, Double>() {
			@Override
			public Double apply(RandomVariable randomVariable) {
				return randomVariable.sub(lnsvqdModel.zeta.apply(randomVariable).mult(deltaT)).sub(c).doubleValue();
			}
		};

		double[] evalPoints = LNSVQDUtils.createTimeGrid(-5, 10, 100);
		for(double value : evalPoints) {
			System.out.println(value + "\t" + rootFunction.apply(lnsvqdModel.getRandomVariableForConstant(value)));
		}
	}

	@Test
	public void printZetaFunction() {
		double[] evalPoints = LNSVQDUtils.createTimeGrid(100, 200, 100);
		for(double value : evalPoints) {
			RandomVariable randomVariable = new Scalar(value);
			RandomVariable randomVariable2 = randomVariable.mult(-1).exp().mult(lnsvqdModel.getKappa1() * lnsvqdModel.getTheta());
			RandomVariable randomVariableZ = randomVariable2.sub(randomVariable2.mult(-1).exp().mult(lnsvqdModel.getKappa2()));
			RandomVariable randomVariable3 = randomVariableZ.add(-lnsvqdModel.getKappa1() + lnsvqdModel.getKappa2() * lnsvqdModel.getTheta()
					- 0.5 * lnsvqdModel.getTotalInstVar());
			double result = randomVariable3.doubleValue();
			System.out.println(value + "\t" + result);
		}
	}

	@Test
	public void checkNewtonsMethod() {
		double root = 5.3804;
		Function<Double, Double> function = new Function<Double, Double>() {
			@Override
			public Double apply(Double aDouble) {
				return 5 * Math.pow(aDouble - root, 2);
			}
		};

		Function<Double, Double> functionDerivatie = new Function<Double, Double>() {
			@Override
			public Double apply(Double aDouble) {
				return 5 * 2 * (aDouble - root);
			}
		};

		NewtonsMethod newtonsMethod = new NewtonsMethod(0);
		for(int k = 0; k < 100; k++) {
			double value = function.apply(newtonsMethod.getNextPoint());
			double derivative = functionDerivatie.apply(newtonsMethod.getNextPoint());
			newtonsMethod.setValueAndDerivative(value, derivative);
		}
		System.out.println("Newton's method root: " + newtonsMethod.getBestPoint());
		Assert.assertEquals(root, newtonsMethod.getBestPoint(), epsilon);
	}

}


