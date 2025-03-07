package net.finmath.equities;

import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDModelCalibrator;
import net.finmath.equities.models.VolatilityPointsSurface;
import org.junit.BeforeClass;
import org.junit.Test;

public class LNSVQDModelCalibratorTest extends TestsSetupForLNSVQD{

	@Test
	public void calibrateTest() throws Exception {
		setDAXHestonSetupCALIB();

		/**
		 * 1. Calibrate and get cvalibrated paramerters
		 */
		double[] calibratedParameters;
		int[] indicesCalibratedParams = {
				0
				/*, 1
				, 2*/
				, 3
				, 4
				, 5
		};
		calibratedParameters = LNSVQDModelCalibrator.calibrate(selectedParamsLNSVQD, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);
		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null);
		impliedVolSurface.printVolSurfaceForOutput();
	}

}