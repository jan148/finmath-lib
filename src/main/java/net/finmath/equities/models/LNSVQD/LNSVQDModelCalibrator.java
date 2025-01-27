package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import org.apache.commons.math3.util.Pair;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

/**
 * The exponential-affine approximation to the mgf can be precalculated; For t0 < t1, we would calculate A_t0 twice,
 * once for A_t0 and once for A_t0.
 */

public class LNSVQDModelCalibrator {

	// TODO
	public static double[] calibrate(final double[] initialVolatilityParameters,
	                                 int[] parameterIndices,
	                                 LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer,
	                                 VolatilityPointsSurface volatilitySurface) throws SolverException {
		LocalDate today = volatilitySurface.getToday();

		// Create list of parameters that should be calibrated
		double[] initialVolatilityParametersToCalibrate = new double[parameterIndices.length];
		for(int j = 0; j < initialVolatilityParametersToCalibrate.length; j++) {
			initialVolatilityParametersToCalibrate[j] = initialVolatilityParameters[parameterIndices[j]];
		}

		double[] calibratedParameters = new double[initialVolatilityParameters.length];

		final double initalStockValue = lnsvqdModelAnalyticalPricer.getSpot0();
		final double riskFreeRate = lnsvqdModelAnalyticalPricer.getRiskFreeRate(1); // Check

		final double[] targetValues = volatilitySurface.getVolatilityPoints().stream()
				.mapToDouble(VolatilityPoint::getVolatility)
				.toArray();
		final List<Pair<Double, Double>> strikeMaturityPairs = volatilitySurface.getVolatilityPoints().stream()
				.map(volaPoint -> {
					double ttm = volatilitySurface.getDayCountConvention().getDaycountFraction(today, volaPoint.getDate());
					return new Pair<>(ttm, volaPoint.getStrike());
				})
				.collect(Collectors.toList());

		// Optimization algorithm parameters
		final int maxIteration = 60;
		final int numberOfThreads = 4;

		/**
		 * For concurrency: Create executor service
		 */
		ExecutorService executorService = new ForkJoinPool(numberOfThreads) ; //

		/**
		 * Create optimizer
		 */
		LevenbergMarquardt levenbergMarquardt = new LevenbergMarquardt(initialVolatilityParametersToCalibrate, targetValues, maxIteration, executorService) {
			@Override
			public void setValues(double[] parameters, double[] values) {
				double[] paramVector = initialVolatilityParameters.clone();
				for(int j = 0; j < parameterIndices.length; j++) {
					int index = parameterIndices[j];
					paramVector[index] = parameters[j];
				}
				lnsvqdModelAnalyticalPricer.setVolatilityParameters(paramVector);

				try {
					double[] impliedVols = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilitySurface).getVolatilityPoints()
							.stream()
							.mapToDouble(VolatilityPoint::getVolatility)
							.toArray();
					System.arraycopy(impliedVols, 0, values, 0, impliedVols.length);
				} catch(Exception e) {
					throw new RuntimeException(e);
				}
			}
		};

		// Set lambda
		levenbergMarquardt.setLambda(0.1);

		// Print some information
		System.out.println("Calibration started");

		// Solve
		levenbergMarquardt.run();

		// Output
		calibratedParameters = levenbergMarquardt.getBestFitParameters();
		return calibratedParameters;
	}

}
