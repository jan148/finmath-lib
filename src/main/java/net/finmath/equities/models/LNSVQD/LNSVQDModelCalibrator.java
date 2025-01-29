package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import net.finmath.time.daycount.DayCountConvention;
import org.apache.commons.math3.fitting.leastsquares.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.lang.Math.E;

/**
 * The exponential-affine approximation to the mgf can be precalculated; For t0 < t1, we would calculate A_t0 twice,
 * once for A_t0 and once for A_t0.
 */

public class LNSVQDModelCalibrator {

	static double shiftSize = 1E-8;

	public static double[] calibrate(final double[] initialVolatilityParameters,
	                                 int[] parameterIndices,
	                                 LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer,
	                                 VolatilityPointsSurface volatilitySurface) throws Exception {
		System.out.println("-------------------------------------");
		final double[] initX = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilitySurface).getVolatilityPoints()
				.stream()
				.mapToDouble(VolatilityPoint::getVolatility)
				.toArray();
		final double[] targetValues = volatilitySurface.getVolatilityPoints().stream()
				.mapToDouble(VolatilityPoint::getVolatility)
				.toArray();
		double initialCost = IntStream.range(0, initX.length)
				.mapToDouble(i -> Math.pow(initX[i] - targetValues[i], 2) * (1. / targetValues.length))
				.sum();
		System.out.println("Calibration started. Initial cost: " + initialCost);
		System.out.print("Initial params: ");
		LNSVQDUtils.printArray(initialVolatilityParameters);

		// Create list of parameters that should be calibrated
		double[] initialVolatilityParametersToCalibrate = new double[parameterIndices.length];
		for(int j = 0; j < initialVolatilityParametersToCalibrate.length; j++) {
			initialVolatilityParametersToCalibrate[j] = initialVolatilityParameters[parameterIndices[j]];
		}

		// Optimization algorithm parameters
		final int maxIteration = 50;
		final int numberOfThreads = 4;

		double[] calibratedParameters = initialVolatilityParameters.clone();
		double[] outputParamsOptimizer = initialVolatilityParameters.clone();

		/**
		 * For concurrency: Create executor service
		 */
		ExecutorService executorService = new ForkJoinPool(numberOfThreads);

		/**
		 * Create optimizer / finmath
		 */
		/*LevenbergMarquardt levenbergMarquardt = new LevenbergMarquardt(
				initialVolatilityParametersToCalibrate
				, targetValues
				, maxIteration
				, executorService) {
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
		levenbergMarquardt.setLambda(0.10);

		// Print some information
		System.out.println("Calibration started");

		// Solve
		levenbergMarquardt.run();

		// Output
		outputParamsOptimizer = levenbergMarquardt.getBestFitParameters();*/

		/**
		 * Create optimizer / Apache
		 */
		// ---------------------------
		LevenbergMarquardtOptimizer levenbergMarquardtOptimizer = new LevenbergMarquardtOptimizer();
		MultivariateJacobianFunction model = params -> {
			double[] paramsFull = initialVolatilityParameters.clone();
			for(int i = 0; i < params.getDimension(); i++) {
				paramsFull[parameterIndices[i]] = params.getEntry(i);
			}
			lnsvqdModelAnalyticalPricer.setVolatilityParameters(paramsFull);
			double[] impliedVols;
			try {
				impliedVols = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilitySurface).getVolatilityPoints()
						.stream()
						.mapToDouble(VolatilityPoint::getVolatility)
						.toArray();
			} catch(Exception e) {
				throw new RuntimeException(e);
			}

			RealVector value = new ArrayRealVector(targetValues.length);
			RealMatrix jacobian = new Array2DRowRealMatrix(targetValues.length, params.getDimension());

			double[][] forwardShiftedValues = new double[targetValues.length][params.getDimension()];
			for(int j = 0; j < params.getDimension(); j++) {
				double[] paramsFullShifted = paramsFull.clone();
				paramsFullShifted[parameterIndices[j]] = paramsFullShifted[parameterIndices[j]] + shiftSize;
				lnsvqdModelAnalyticalPricer.setVolatilityParameters(paramsFullShifted);
				try {
					forwardShiftedValues[j] = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilitySurface).getVolatilityPoints()
							.stream()
							.mapToDouble(VolatilityPoint::getVolatility)
							.toArray();
				} catch(Exception e) {
					throw new RuntimeException(e);
				}
				for(int i = 0; i < forwardShiftedValues[j].length; i++) {
					forwardShiftedValues[j][i] = (forwardShiftedValues[j][i] - impliedVols[i]) / shiftSize; // Divide each element by shiftSize
				}
			}

			for(int i = 0; i < targetValues.length; i++) {
				value.setEntry(i, impliedVols[i]);
				for(int j = 0; j < params.getDimension(); j++) {
					jacobian.setEntry(i, j, forwardShiftedValues[j][i]);
				}
			}

			return new Pair<>(value, jacobian);
		};
		LeastSquaresProblem leastSquaresProblem = new LeastSquaresBuilder()
				.start(initialVolatilityParametersToCalibrate)
				.model(model)
				.target(targetValues)
				.lazyEvaluation(false)
				.maxEvaluations(50)
				.maxIterations(maxIteration)
				.build();
		// Run the optimizer
		outputParamsOptimizer = levenbergMarquardtOptimizer.optimize(leastSquaresProblem).getPoint().toArray();

		/**
		 * Get summary
		 */
		// Retrieve the optimal parameters
		for(int i = 0; i < parameterIndices.length; i++) {
			calibratedParameters[parameterIndices[i]] = outputParamsOptimizer[i];
		}

		// Retrieve the optimal value (cost function value)
		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		double[] endX = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilitySurface).getVolatilityPoints()
				.stream()
				.mapToDouble(VolatilityPoint::getVolatility)
				.toArray();
		double finalCost = IntStream.range(0, initX.length)
				.mapToDouble(i -> Math.pow(endX[i] - targetValues[i], 2) * (1. / targetValues.length))
				.sum();
		System.out.println("\nCalibration ended. Final cost: " + finalCost);
		System.out.print("Final params:");
		LNSVQDUtils.printArray(calibratedParameters);
		System.out.println("\n-------------------------------------");
		// ---------------------------

		return calibratedParameters;
	}
}
