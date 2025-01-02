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
import java.util.function.BiFunction;

public class LNSVQDModelCalibrator {

	// TODO
	public static double[] calibrate(final double[] initialVolatilityParameters, LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer,
	                             DynamicVolatilitySurface volatilitySurface) throws SolverException {
		double[] calibratedParameters = new double[initialVolatilityParameters.length];

		final double initalStockValue = lnsvqdModelAnalyticalPricer.getSpot0();
		final double riskFreeRate = lnsvqdModelAnalyticalPricer.getRiskFreeRate();

		final double[] targetValues = volatilitySurface.getPrices(initalStockValue, riskFreeRate)
				.stream().mapToDouble(Double::doubleValue).toArray();
		final int maxIteration = 5;
		final int numberOfThreads = 1;

		final LocalDate calibrationDay = volatilitySurface.getToday();

		LevenbergMarquardt levenbergMarquardt = new LevenbergMarquardt(initialVolatilityParameters, targetValues, maxIteration, numberOfThreads) {
			@Override
			public void setValues(double[] parameters, double[] values) throws SolverException {
				// Set parameters
				lnsvqdModelAnalyticalPricer.setVolatilityParameters(parameters);

				for(int p = 0; p < volatilitySurface.getNumberOfVolatilityPoints(); p++) {
					// Get ttm
					LocalDate date = volatilitySurface.getVolatilityDates().get(p);
					long period = ChronoUnit.DAYS.between(calibrationDay, date);
					double ttm = period / 365; // TODO: Get the correct denominator

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

		// Print some information
		System.out.println("Calibration started");

		// Solve
		levenbergMarquardt.run();

		// Output
		calibratedParameters = levenbergMarquardt.getBestFitParameters();
		return calibratedParameters;
	}

}
