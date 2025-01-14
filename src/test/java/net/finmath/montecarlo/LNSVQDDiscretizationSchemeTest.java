package net.finmath.montecarlo;

import net.finmath.equities.models.LNSVQD.LNSVQDModel;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
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
	int numberOfPaths = 100;
	int seed = 4235;

	/**
	 * Model params
	 */
	// Right params: sigma0=0.8327, theta=1.0139, kappa1=4.8606, kappa2=4.7938, beta=0.1985, volvol=2.3690
	private final double spot0 = 1;
	private final double sigma0 = 0.8327; //0.41;
	// Value as in paper
	private final double kappa1 = 4.8606; // 4.8606; //2.21;
	// Value as in paper
	private final double kappa2 = 4.7938; // 4.7938; //2.18;
	private final double theta = 1.0139; // 1.0139; //0.38;
	private final double beta = 0.1985; // 0.1985; //0.5;
	private final double epsilon = 2.3690; //2.3690; //3.06;


	/**
	 * Models
	 */
	LNSVQDModel lnsvqdModel = new LNSVQDModel(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);

	/**
	 * Option params
	 */
	double strike = 1;
	double maturity = 1.2; // 1.2

	/**
	 * Market observables
	 */
	double riskFreeRate = lnsvqdModel.getRiskFreeRate();
	double discountFactor = Math.exp(-riskFreeRate * maturity);
	double convenienceFcator = 0;

	/**
	 * Time discretization
	 */
	double[] timeGrid = LNSVQDUtils.createTimeGrid(0, maturity, 100);
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
		int componentIndex = 0;

		// Create 2D-Brownian motion and discretization scheme
		BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
		LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdModel, brownianMotion);

		for(int i = 0; i < numberOfPaths; i++) {
			System.out.print("\t" + "Realization " + i);
		}
		for(int k = 0; k < timeDiscretization.getNumberOfTimes(); k++) {
			System.out.print("\n" + timeDiscretization.getTime(k) + "\t");
			for(int j = 0; j < numberOfPaths; j++) {
				System.out.print(lnsvqdDiscretizationScheme.getProcessValue(k, componentIndex).get(j) + "\t");
			}
		}
		System.out.print("\n\n");
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


