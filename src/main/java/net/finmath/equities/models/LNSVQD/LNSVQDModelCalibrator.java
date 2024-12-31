package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.DynamicVolatilitySurface;
import net.finmath.equities.models.SviVolatilitySmile;
import net.finmath.equities.models.SviVolatilitySurface;
import net.finmath.equities.models.VolatilitySurface;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import org.apache.commons.math3.complex.Complex;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;

public class LNSVQDModelCalibrator {
	boolean isCalibrated = false;
	/**
	 * Numerical parameters
	 */
	/*private final int numStepsForODEIntegration = 100;
	private final int numStepsForInfiniteIntegral = 1000;
	private final double upperBoundForInfiniteIntegral = numStepsForInfiniteIntegral / 10;*/

	/**
	 * Inital parameters
	 */
	final double[] initialVolatilityParameters;

	public LNSVQDModelCalibrator(double[] initialVolatilityParameters) {
		this.initialVolatilityParameters = initialVolatilityParameters;
	}

	// TODO
	public void calibrate(LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer, DynamicVolatilitySurface volatilitySurface) {
		double loss = 0;

		final double[] initialParameters = initialVolatilityParameters;
		final double[] targetValues = volatilitySurface.getVolatilityPoints().stream()
				.mapToDouble(volatilityPoint -> volatilityPoint.getVolatility()).toArray();
		final int maxIteration = 500;
		final int numberOfThreads = 1;

		LevenbergMarquardt levenbergMarquardt = new LevenbergMarquardt(initialVolatilityParameters, targetValues, maxIteration, numberOfThreads) {
			@Override
			public void setValues(double[] parameters, double[] values) throws SolverException {
				// Set parameterss
				lnsvqdModelAnalyticalPricer.setVolatilityParameters(parameters);

				for(int p = 0; p < volatilitySurface.getNumberOfVolatilityPoints(); p++) {
					// Fetch point from the surface
					LocalDate date = volatilitySurface.getVolatilityDates().get(p);
					double strike = volatilitySurface.getStrikes().get(p);
					double empiricalVolatility = volatilitySurface.getVolatility(strike, date, null);
					// double empiricalBSCallPrice = func(empiricalVolatility, ...);

					// Get ttm
					double ttm = 0; // TODO
					double discountFactor = Math.exp(-lnsvqdModelAnalyticalPricer.getRiskFreeRate() * ttm);
					double conveniceFacor = 0; // TODO: Change later

					// Calculate corresponding call option price
					double lnsvqdCallPrice = lnsvqdModelAnalyticalPricer.getCallPrice(strike, ttm, discountFactor, conveniceFacor);

					// Value
					values[p] = lnsvqdCallPrice - empiricalBSCallPrice;
				}
			}
		};

		for(int p = 0; p < volatilitySurface.getNumberOfVolatilityPoints(); p++) {
			// Fetch point from the surface
			LocalDate date = volatilitySurface.getVolatilityDates().get(p);
			double strike = volatilitySurface.getStrikes().get(p);
			double empiricalVolatility = volatilitySurface.getVolatility(strike, date, null);

			// Calculate
		}



		// int numberDates = volatilitySurface.length;

	}

}
