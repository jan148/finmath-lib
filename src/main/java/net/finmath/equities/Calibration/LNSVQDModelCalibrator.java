package net.finmath.equities.Calibration;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.equities.pricer.LNSVQDModelAnalyticalPricer;
import net.finmath.marketdata.model.curves.DiscountCurve;
import org.apache.commons.math3.fitting.leastsquares.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.Pair;

import java.time.LocalDate;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * This class calibrates the LNSVQD to a given IVS.
 *
 * @author Jan Berger
 */

public class LNSVQDModelCalibrator {
	private static double shiftSize = 1E-8;

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

	public static double[] calibrate(final double[] initialVolatilityParameters,
	                                 int[] parameterIndices,
	                                 final LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer,
	                                 VolatilityPointsSurface volatilitySurface) throws Exception {
		double spot0 = lnsvqdModelAnalyticalPricer.getSpot0();
		LocalDate valuationDate = lnsvqdModelAnalyticalPricer.getSpotDate();
		YieldCurve discountCurve = lnsvqdModelAnalyticalPricer.getDiscountCurve();
		EquityForwardStructure equityForwardStructure = lnsvqdModelAnalyticalPricer.getEquityForwardStructure();
		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricerMain = lnsvqdModelAnalyticalPricer.copyWithNewParameters(
				spot0, initialVolatilityParameters[0], initialVolatilityParameters[1], initialVolatilityParameters[2], initialVolatilityParameters[3]
				, initialVolatilityParameters[4], initialVolatilityParameters[5], 0, valuationDate, discountCurve, equityForwardStructure);

		System.out.println("-------------------------------------");
		final double[] initX = lnsvqdModelAnalyticalPricerMain.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
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

			LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricerNew = lnsvqdModelAnalyticalPricer.copyWithNewParameters(
					spot0, paramsFull[0], paramsFull[1], paramsFull[2], paramsFull[3], paramsFull[4], paramsFull[5], 0
					, valuationDate, discountCurve, equityForwardStructure);
			double[] impliedVols;
			try {
				impliedVols = lnsvqdModelAnalyticalPricerNew.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
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
				LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricerShifted = lnsvqdModelAnalyticalPricer.copyWithNewParameters(
						spot0, paramsFullShifted[0], paramsFullShifted[1], paramsFullShifted[2], paramsFullShifted[3], paramsFullShifted[4], paramsFullShifted[5], 0
						, valuationDate, discountCurve, equityForwardStructure);
				try {
					forwardShiftedValues[j] = lnsvqdModelAnalyticalPricerShifted.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
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

		// Retrieve the optimal parameters
		for(int i = 0; i < parameterIndices.length; i++) {
			calibratedParameters[parameterIndices[i]] = outputParamsOptimizer[i];
		}

		// Summary
		// Retrieve the optimal value (cost function value)
		lnsvqdModelAnalyticalPricerMain = lnsvqdModelAnalyticalPricer.copyWithNewParameters(
				spot0, calibratedParameters[0], calibratedParameters[1], calibratedParameters[2]
				, calibratedParameters[3], calibratedParameters[4], calibratedParameters[5], 0
				, valuationDate, discountCurve, equityForwardStructure);
		double[] endX = lnsvqdModelAnalyticalPricerMain.getImpliedVolSurfaceFromVolSurface(volatilitySurface, null).getVolatilityPoints()
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
