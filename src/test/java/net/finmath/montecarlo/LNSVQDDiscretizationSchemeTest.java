package net.finmath.montecarlo;

import net.finmath.equities.models.LNSVQDModel;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.LNSVQDDiscretizationScheme;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import org.junit.Test;

import static org.junit.jupiter.api.Assertions.*;

public class LNSVQDDiscretizationSchemeTest {
	/**
	 * Simulation parameters
	 */
	int numberOfPaths = 1000;
	int seed = 3107;

	/**
	 * Model params
	 */
	private final double spot0 = 1;
	private final double sigma0 = 0.5;
	private final double kappa1 = 0;
	private final double kappa2 = 0;
	private final double theta = 0;
	private final double beta = 0;
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
	double[] timeGrid = LNSVQDUtils.createTimeGrid(0, 10, 20);
	TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(timeGrid);

	/**
	 * Time discretization
	 */
	RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

	/**
	 * Tolerance level
	 */
	private final double delta =  10E-3;

	@Test
	public void doPrecalculateProcess() throws CalculationException {
		// Print time discretization
		// for(double t : timeDiscretization.getAsDoubleArray()) {System.out.println(t);}

		// Create 2D-Brownian motion
		BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
		RandomVariable brownianIncrementAt0 = brownianMotion.getBrownianIncrement(0, 0);
		/*for(int j = 0; j < numberOfPaths; j++) {
			System.out.println(brownianIncrementAt0.get(j));
		}*/

		LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdModel, brownianMotion);
		System.out.println(lnsvqdDiscretizationScheme.getProcessValue(0, 1));

	}



}