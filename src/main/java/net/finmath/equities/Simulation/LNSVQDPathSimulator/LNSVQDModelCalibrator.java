package net.finmath.equities.Simulation.LNSVQDPathSimulator;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.equities.pricer.LNSVQDModelAnalyticalPricer;
import org.apache.commons.math3.fitting.leastsquares.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

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
		final double[] initX = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
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
		System.out.println("-------------------------------------");

		// Create list of parameters that should be calibrated
		double[] initialVolatilityParametersToCalibrate = new double[parameterIndices.length];
		for(int j = 0; j < initialVolatilityParametersToCalibrate.length; j++) {
			initialVolatilityParametersToCalibrate[j] = initialVolatilityParameters[parameterIndices[j]];
		}

		// Optimization algorithm parameters
		final int maxIteration = 100;
		final int maxEvaluations = 200;
		AtomicInteger iterationCount = new AtomicInteger(0);

		double[] calibratedParameters = initialVolatilityParameters.clone();
		double[] outputParamsOptimizer = initialVolatilityParameters.clone();

		/**
		 * Create optimizer / Apache
		 */
		// ---------------------------
		LevenbergMarquardtOptimizer levenbergMarquardtOptimizer = new LevenbergMarquardtOptimizer()
				.withCostRelativeTolerance(1.0e-6)
				.withParameterRelativeTolerance(1.0e-6)
				.withOrthoTolerance(1.0e-6);

		MultivariateJacobianFunction model = params -> {
			int iter = iterationCount.incrementAndGet(); // Atomic and thread-safe
			System.out.println("Iteration: " + iter);

			double[] paramsFull = initialVolatilityParameters.clone();
			for(int i = 0; i < params.getDimension(); i++) {
				paramsFull[parameterIndices[i]] = params.getEntry(i);
			}

			System.out.print("Current params: ");
			LNSVQDUtils.printArray(Arrays.stream(paramsFull).toArray());

			lnsvqdModelAnalyticalPricer.setVolatilityParameters(paramsFull);
			double[] impliedVols;
			try {
				impliedVols = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
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
					forwardShiftedValues[j] = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
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

			double mse = IntStream.range(0, impliedVols.length)
					.mapToDouble(i -> Math.pow(impliedVols[i] - targetValues[i], 2) * (1. / impliedVols.length))
					.sum();
			System.out.println("MSE before iteration " + iter + ": " + mse);

			return new Pair<>(value, jacobian);
		};
		LeastSquaresProblem leastSquaresProblem = new LeastSquaresBuilder()
				.start(initialVolatilityParametersToCalibrate)
				.model(model)
				.target(targetValues)
				.lazyEvaluation(false)
				.maxEvaluations(maxEvaluations)
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
		double[] endX = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
				.stream()
				.mapToDouble(VolatilityPoint::getVolatility)
				.toArray();
		double finalCost = IntStream.range(0, initX.length)
				.mapToDouble(i -> Math.pow(endX[i] - targetValues[i], 2) * (1. / targetValues.length))
				.sum();
		System.out.println("\nCalibration ended. Final cost: " + finalCost);
		System.out.print("Final params: ");
		LNSVQDUtils.printArray(calibratedParameters);
		System.out.println("\n-------------------------------------");
		// ---------------------------

		return calibratedParameters;
	}
}
