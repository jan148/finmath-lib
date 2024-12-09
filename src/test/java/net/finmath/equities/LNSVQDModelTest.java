package net.finmath.equities;

import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.LNSVQDModel;
import net.finmath.functions.AnalyticFormulas;
import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * This test compares the LNSVQD call option price with the BS call option price
 * TODO: Extend to Heston call option price
 */
class LNSVQDModelTest {
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
	double maturity = 1;

	/**
	 * Market observables
	 */
	double riskFreeRate = 0.05;
	double discountFactor = Math.exp(-riskFreeRate * maturity);
	double convenienceFcator = 0;

	@Test
	void calculateExponentialAffineApproximation() {
		double y = 1;
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
		Complex exponentialAffineApproximationOdeValue = lnsvqdModel.calculateExponentialAffineApproximation(maturity, charFuncArgs);
		Complex exponentialAffineApproximationAnalyticalValue = (charFuncArgs[0].multiply(-lnsvqdModel.getX0())
				.add(charFuncArgs[0].pow(2).add(charFuncArgs[0]).subtract(charFuncArgs[1]).multiply(0.5 * maturity).multiply(Math.pow(lnsvqdModel.getY0(), 2))))
				.exp();
		System.out.println("Analytical exponential-affine approximation value at " + maturity + ": " + exponentialAffineApproximationAnalyticalValue);
		System.out.println("ODE-based exponential-affine approximation value at " + maturity + ": " + exponentialAffineApproximationOdeValue);
	}

	@Test
	void getCallPrice() {
		//TODO: Implement forward structure

		// Get option values
		double bsOptionValue = AnalyticFormulas.blackScholesOptionValue(spot0, riskFreeRate, sigma0, maturity, strike, true);
		double lnsvqdOptionValue = lnsvqdModel.getCallPrice(strike, maturity, discountFactor, convenienceFcator);

		// Print
		System.out.println("Call oprion price BS: \t" + bsOptionValue);
		System.out.println("Call oprion price LNSVQD: \t" + lnsvqdOptionValue);
	}


}
