package net.finmath.equities;

import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDModelCalibrator;
import net.finmath.equities.models.VolatilityPointsSurface;
import org.apache.commons.math3.util.Pair;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.ArrayList;

public class LNSVQDModelCalibratorTest extends TestsSetupForLNSVQD{

	@Test
	public void calibrateTest() throws Exception {
		setDAXHestonFebruarySetupSIM();

		/**
		 * 1. Calibrate and get cvalibrated paramerters
		 */
		double[] calibratedParameters;
		int[] indicesCalibratedParams = {
				0
				/*, 1*/
				/*, 2*/
				, 3
				, 4
				, 5
		};
		calibratedParameters = LNSVQDModelCalibrator.calibrate(selectedParamsToCalibrate, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);
		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null);
		impliedVolSurface.printVolSurfaceForOutput();
	}

}