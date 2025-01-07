package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.DynamicVolatilitySurface;
import net.finmath.equities.models.SviVolatilitySmile;
import net.finmath.equities.models.SviVolatilitySurface;
import net.finmath.equities.models.VolatilitySurface;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import org.apache.commons.math3.complex.Complex;

import java.time.LocalDate;
import java.time.Period;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.function.BiFunction;

/**
 * The exponential-affine approximation to the mgf can be precalculated; For t0 < t1, we would calculate A_t0 twice,
 * once for A_t0 and once for A_t0.
 */

public class LNSVQDModelCalibrator {

	// TODO
	public static double[] calibrate(final double[] initialVolatilityParameters, LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer,
	                             DynamicVolatilitySurface volatilitySurface, int[] parameterIndices) throws SolverException {
		// Create list of paramters that should be calibrated
		double[] initialVolatilityParametersToCalibrate = new double[parameterIndices.length];
		for(int j = 0; j < initialVolatilityParametersToCalibrate.length; j++) {
			initialVolatilityParametersToCalibrate[j] = initialVolatilityParameters[parameterIndices[j]];
		}

		double[] calibratedParameters = new double[initialVolatilityParameters.length];

		final double initalStockValue = lnsvqdModelAnalyticalPricer.getSpot0();
		final double riskFreeRate = lnsvqdModelAnalyticalPricer.getRiskFreeRate();

		final double[] targetValues = volatilitySurface
				.getPrices(initalStockValue, riskFreeRate, volatilitySurface.getDayCountConvention())
				.stream().mapToDouble(Double::doubleValue).toArray();
		final int maxIteration = 60;
		final int numberOfThreads = 1;

		final LocalDate calibrationDay = volatilitySurface.getToday();

		/**
		 * For concurrency: Create executor service
		 */
		ExecutorService executorService = null ; // new ForkJoinPool(numberOfThreads)

		/**
		 * Create optimizer
		 */
		LevenbergMarquardt levenbergMarquardt = new LevenbergMarquardt(initialVolatilityParametersToCalibrate, targetValues, maxIteration, executorService) {
			@Override
			public void setValues(double[] parameters, double[] values) throws SolverException {
				double[] paramVector = initialVolatilityParameters.clone();
				for(int j = 0; j < parameterIndices.length; j++) {
					int index = parameterIndices[j];
					paramVector[index] = parameters[j];
				}

				// Set parameters
				lnsvqdModelAnalyticalPricer.setVolatilityParameters(paramVector);

				for(int p = 0; p < volatilitySurface.getNumberOfVolatilityPoints(); p++) {
					// Get ttm
					LocalDate date = volatilitySurface.getVolatilityDates().get(p);
					long period = ChronoUnit.DAYS.between(calibrationDay, date);
					double ttm = period / 365.; // TODO: Get the correct denominator

					double discountFactor = Math.exp(-lnsvqdModelAnalyticalPricer.getRiskFreeRate() * ttm);
					double conveniceFacor = 0; // TODO: Change later

					// Fetch point from the surface
					double strike = volatilitySurface.getStrikes().get(p);

					// Calculate corresponding call option price
					double lnsvqdCallPrice = lnsvqdModelAnalyticalPricer.getCallPrice(strike, ttm, discountFactor, conveniceFacor);

					// Value
					values[p] = lnsvqdCallPrice;
				}
			}
		};

		// Set lambda
		levenbergMarquardt.setLambda(0.001);

		// Print some information
		System.out.println("Calibration started");

		// Solve
		levenbergMarquardt.run();

		// Output
		calibratedParameters = levenbergMarquardt.getBestFitParameters();
		return calibratedParameters;
	}

}
