package net.finmath.equities;

import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.LNSVQDModel;
import net.finmath.functions.AnalyticFormulas;
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
	double riskFreeRate = 0.1;
	double discountFactor = Math.exp(-riskFreeRate * maturity);
	double convenienceFcator = 0;

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