package net.finmath.montecarlo.assetderivativevaluation.models;

import java.util.Map;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloAssetModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

/**
 *
 * TODO Add UnitTest for correlation - using ExchangeOption.
 *
 * @author Christian Fries
 */
public class MultiAssetBlackScholesModelTest {

	@Test
	public void testModelWithFactorLoadings() throws CalculationException {

		double[] initialValue = new double[] { 100.0, 110 };
		double riskFreeRate = 0.05;
		double[][] factorLoadings = new double[][] {
			{ 0.25, 0.10 },
			{ 0.00, 0.20 }

		};

		MultiAssetBlackScholesModel model = new MultiAssetBlackScholesModel(initialValue, riskFreeRate, factorLoadings);

		int numberOfPaths = 2000000;
		int seed = 3216;
		MonteCarloAssetModel assetModel = getMonteCarloAssetModel(model, numberOfPaths, seed);

		/*
		 * Test the two individual European options
		 */
		double maturity = 5.0;
		double strike1 = 120;
		double strike2 = 130;
		EuropeanOption europeanOption1 = new EuropeanOption(maturity, strike1, 0);
		EuropeanOption europeanOption2 = new EuropeanOption(maturity, strike2, 1);

		double volatility1 = Math.sqrt(factorLoadings[0][0]*factorLoadings[0][0] + factorLoadings[0][1]*factorLoadings[0][1]);

		assertEqualsValuationEuropeanProduct(assetModel, europeanOption1, initialValue[0], riskFreeRate, volatility1);

		double volatility2 = Math.sqrt(factorLoadings[1][0]*factorLoadings[1][0] + factorLoadings[1][1]*factorLoadings[1][1]);
		assertEqualsValuationEuropeanProduct(assetModel, europeanOption2, initialValue[1], riskFreeRate, volatility2);
	}

	@Test
	void testModelWithVolatilityAndCorrelation() throws CalculationException {

		double[] initialValue = new double[] { 100.0, 110 };
		double riskFreeRate = 0.05;
		double[] volatilities = new double[] { 0.30, 0.20 };
		double[][] correlations = new double[][] {
			{ 1.0, 0.0 },
			{ 0.0, 1.0 }

		};

		MultiAssetBlackScholesModel model = new MultiAssetBlackScholesModel(initialValue, riskFreeRate, volatilities, correlations);

		int numberOfPaths = 2000000;
		int seed = 3216;
		MonteCarloAssetModel assetModel = getMonteCarloAssetModel(model, numberOfPaths, seed);

		/*
		 * Test the two individual European options
		 */
		double maturity = 5.0;
		double strike1 = 120;
		double strike2 = 130;
		EuropeanOption europeanOption1 = new EuropeanOption(maturity, strike1, 0);
		EuropeanOption europeanOption2 = new EuropeanOption(maturity, strike2, 1);

		assertEqualsValuationEuropeanProduct(assetModel, europeanOption1, initialValue[0], riskFreeRate, volatilities[0]);

		assertEqualsValuationEuropeanProduct(assetModel, europeanOption2, initialValue[1], riskFreeRate, volatilities[1]);
	}

	private MonteCarloAssetModel getMonteCarloAssetModel(MultiAssetBlackScholesModel model, int numberOfPaths, int seed) {
		int numberOfTimeSteps = 30;
		double deltaT = 0.2;

		int numberOfFactors = model.getNumberOfFactors();

		// Create a time discretization
		final TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0 /* initial */, numberOfTimeSteps, deltaT);

		// Create Brownian motion (and random number generator)
		final BrownianMotion brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, numberOfFactors , numberOfPaths, seed);

		MonteCarloProcess process = new EulerSchemeFromProcessModel(model, brownianMotion);

		return new MonteCarloAssetModel(process);
	}

	@Test
	public void testModelCloneWithModifiedSeed() throws CalculationException {

		double[] initialValue = new double[] { 100.0, 110 };
		double riskFreeRate = 0.05;
		double[][] factorLoadings = new double[][] {
			{ 0.25, 0.10 },
			{ 0.00, 0.20 }

		};

		MultiAssetBlackScholesModel model = new MultiAssetBlackScholesModel(initialValue, riskFreeRate, factorLoadings);

		int numberOfPaths = 10000;
		int seed = 3216;
		MonteCarloAssetModel assetModel = getMonteCarloAssetModel(model, numberOfPaths, seed);

		/*
		 * Test the European option with different seeds.
		 */
		double maturity = 5.0;
		double strike = 120;
		EuropeanOption europeanOption = new EuropeanOption(maturity, strike, 0);

		double valueSeed1 = europeanOption.getValue(assetModel);
		double valueSeed2 = europeanOption.getValue(assetModel.getCloneWithModifiedData(Map.of("seed", 3141)));
		double valueSeed3 = europeanOption.getValue(getMonteCarloAssetModel(model, numberOfPaths, 3141));
		
		Assertions.assertNotEquals(valueSeed1, valueSeed2, "Models with different seed.");
		Assertions.assertEquals(valueSeed2, valueSeed3, "Models with same seed, different construction.");
	}

	private static void assertEqualsValuationEuropeanProduct(MonteCarloAssetModel assetModel, EuropeanOption europeanOption, double initialValue, double riskFreeRate, double volatility) throws CalculationException {

		double valueMonteCarlo = europeanOption.getValue(assetModel);
		double valueAnalytic = AnalyticFormulas.blackScholesOptionValue(initialValue, riskFreeRate, volatility, europeanOption.getMaturity(), europeanOption.getStrike());
		Assertions.assertEquals(valueAnalytic, valueMonteCarlo, 5E-2, "valuation");

	}
}
