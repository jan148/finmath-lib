package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.models.DynamicVolatilitySurface;
import net.finmath.equities.models.SviVolatilitySmile;
import net.finmath.equities.models.SviVolatilitySurface;
import net.finmath.equities.models.VolatilitySurface;
import org.apache.commons.math3.complex.Complex;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

public class LNSVQDModelCalibrator {
	/**
	 * Numerical parameters
	 */
	private final int numStepsForODEIntegration = 100;
	private final int numStepsForInfiniteIntegral = 1000;
	private final double upperBoundForInfiniteIntegral = numStepsForInfiniteIntegral / 10;

	public LNSVQDModelCalibrator() {}

	// TODO
	public static void calibrate(LNSVQDModel lnsvqdModel, DynamicVolatilitySurface volatilitySurface) {
		double loss = 0;

		/**
		 * Get LNSVQDModel
		 */


		// int numberDates = volatilitySurface.length;



		for(int s = 0; s < numberDates; s++) {
			// Fetch point from the surface
			LocalDate date = volatilitySurface.getVolatilityDates().get(s);
			double strike = volatilitySurface.getStrikes().get(s);
			double empiricalVolatility = volatilitySurface.getVolatility(strike, date, null);

			// Calculate
		}
	}

}
