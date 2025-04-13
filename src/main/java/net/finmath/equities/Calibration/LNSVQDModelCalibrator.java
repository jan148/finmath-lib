package net.finmath.equities.Calibration;

import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulator;
import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.equities.pricer.LNSVQDModelAnalyticalPricer;
import org.apache.commons.math3.fitting.leastsquares.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.Pair;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * The exponential-affine approximation to the mgf can be precalculated; For t0 < t1, we would calculate A_t0 twice,
 * once for A_t0 and once for A_t0.
 */

public class LNSVQDModelCalibrator {
	static double shiftSize = 1E-8;


	public static double[] getAcfsFromVolPaths(double[][] volaPaths, int numLags, int lagSize) {
		PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation();
		int numberOfPaths = volaPaths[0].length;
		double[] empiricalAcfAtLags = new double[numLags];
		double[][] empiricalAcfAtLagsAndPaths = new double[numLags][numberOfPaths];
		try {
			for(int p = 0; p < volaPaths[0].length; p++) {
				int pathNumber = p;
				double[] volaPath = IntStream.range(0, volaPaths.length)
						.mapToDouble(i -> volaPaths[i][pathNumber]).toArray();
				for(int l = 0; l < numLags; l++) {
					double[] volaPathShifted = Arrays.copyOfRange(volaPath, l * lagSize, volaPath.length);
					double[] volaPathShortened = Arrays.copyOfRange(volaPath, 0, volaPathShifted.length);
					empiricalAcfAtLagsAndPaths[l][p] = pearsonsCorrelation.correlation(volaPathShortened, volaPathShifted);
				}
			}
		} catch(Exception e) {
			throw new RuntimeException(e);
		}
		for(int k = 0; k < numLags; k++) {
			empiricalAcfAtLags[k] = Arrays.stream(empiricalAcfAtLagsAndPaths[k]).average().getAsDouble();
		}
		return empiricalAcfAtLags;
	}

	// Calibrate mean reversion params
	public static double[] calibrateMRParams(final double[] initialVolatilityParameters, LNSVQDPathSimulator lnsvqdPathSimulator,
	                                         List<Pair<Double, Double>> acfAtLags) {
		int lagSize = 10; // 10 days lags
		int numLags = acfAtLags.size();

		System.out.println("------- CALIBRATION OF MEAN-REVERSION PARAMETERS -------");
		// assert(initialVolatilityParameters[0] == 0 && initialVolatilityParameters[4] == 0) : "simga0 & beta need to be set to zero for calibration of mr-params!";
		int[] parameterIndices = new int[]{1, 2};
		lnsvqdPathSimulator.setVolatilityParameters(initialVolatilityParameters);
		lnsvqdPathSimulator.precalculatePaths(3105, false, 0, null, Boolean.TRUE);
		double[][] volPaths = lnsvqdPathSimulator.path[1];
		double[] initX = getAcfsFromVolPaths(volPaths, numLags, lagSize);
		final double[] targetValues = acfAtLags.stream()
				.mapToDouble(x -> x.getValue())
				.toArray();
		double initialCost = IntStream.range(0, initX.length)
				.mapToDouble(i -> Math.pow(initX[i] - targetValues[i], 2) * (1. / targetValues.length))
				.sum();
		System.out.println("Calibration started. Initial cost: " + initialCost);
		System.out.print("Initial params: ");
		LNSVQDUtils.printArray(initialVolatilityParameters);
		System.out.println("\n-------------------------------------");

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
			int iter = iterationCount.incrementAndGet();
			System.out.println("Iteration: " + iter);

			double[] paramsFull = initialVolatilityParameters.clone();
			for(int i = 0; i < params.getDimension(); i++) {
				paramsFull[parameterIndices[i]] = params.getEntry(i);
			}

			System.out.print("Current params: ");
			LNSVQDUtils.printArray(Arrays.stream(paramsFull).toArray());

			lnsvqdPathSimulator.setVolatilityParameters(paramsFull);
			lnsvqdPathSimulator.precalculatePaths(3105, false, 0, null, Boolean.TRUE);
			double[][] volPathsCurrentParams = lnsvqdPathSimulator.path[1];
			double[] modelAcf = getAcfsFromVolPaths(volPathsCurrentParams, numLags, lagSize);

			RealVector value = new ArrayRealVector(targetValues.length);
			RealMatrix jacobian = new Array2DRowRealMatrix(targetValues.length, params.getDimension());

			double[][] forwardShiftedValues = new double[params.getDimension()][targetValues.length];
			for(int j = 0; j < params.getDimension(); j++) {
				double[] paramsFullShifted = paramsFull.clone();
				paramsFullShifted[parameterIndices[j]] = paramsFullShifted[parameterIndices[j]] + shiftSize;
				lnsvqdPathSimulator.setVolatilityParameters(paramsFullShifted);
				lnsvqdPathSimulator.precalculatePaths(3105, false, 0, null, Boolean.TRUE);
				double[][] volPathsAtShiftedParams = lnsvqdPathSimulator.path[1];
				double[] modelAcfAtShiftedParams = getAcfsFromVolPaths(volPathsAtShiftedParams, numLags, lagSize);
				for(int i = 0; i < forwardShiftedValues[j].length; i++) {
					forwardShiftedValues[j][i] = (modelAcfAtShiftedParams[i] - modelAcf[i]) / shiftSize; // Divide each element by shiftSize
				}
			}

			for(int i = 0; i < targetValues.length; i++) {
				value.setEntry(i, modelAcf[i]);
				for(int j = 0; j < params.getDimension(); j++) {
					jacobian.setEntry(i, j, forwardShiftedValues[j][i]);
				}
			}

			double mse = IntStream.range(0, modelAcf.length)
					.mapToDouble(i -> Math.pow(modelAcf[i] - targetValues[i], 2) * (1. / modelAcf.length))
					.sum();
			System.out.println("MSE before iteration " + iter + ": " + mse);

			return new Pair<>(value, jacobian);
		};
		ParameterValidator parameterValidator = new ParameterValidator() {
			@Override
			public RealVector validate(RealVector params) {
				RealVector paramsNew = params;
				if(params.getEntry(0) < 0) {
					paramsNew.setEntry(0, Math.abs(params.getEntry(0)));
				}
				if(params.getEntry(0) < 0) {
					paramsNew.setEntry(0, Math.abs(params.getEntry(1)));
				}
				if(params.getEntry(0) > 10) {
					paramsNew.setEntry(0, Math.random() * 5);
				}
				if(params.getEntry(1) > 10) {
					paramsNew.setEntry(1, Math.random() * 5);
				}
				return paramsNew;
			}
		};

		LeastSquaresProblem leastSquaresProblem = new LeastSquaresBuilder()
				.start(initialVolatilityParametersToCalibrate)
				.model(model)
				.target(targetValues)
				.lazyEvaluation(false)
				.maxEvaluations(maxEvaluations)
				.maxIterations(maxIteration)
				.parameterValidator(parameterValidator)
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
		lnsvqdPathSimulator.setVolatilityParameters(calibratedParameters);
		lnsvqdPathSimulator.precalculatePaths(3105, false, 0, null, Boolean.TRUE);
		double[][] volPathsFinal = lnsvqdPathSimulator.path[1];
		double[] endX = getAcfsFromVolPaths(volPathsFinal, numLags, lagSize);
		LNSVQDUtils.printArray(endX);
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

	public static double[] calibrate(final double[] initialVolatilityParameters,
	                                 int[] parameterIndices,
	                                 LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer,
	                                 VolatilityPointsSurface volatilitySurface) throws Exception {
		System.out.println("-------------------------------------");
		lnsvqdModelAnalyticalPricer.setVolatilityParameters(initialVolatilityParameters);
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
		System.out.println("\n-------------------------------------");

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
			int iter = iterationCount.incrementAndGet();
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

			double[][] forwardShiftedValues = new double[params.getDimension()][targetValues.length];
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
