package net.finmath.equities;

import net.finmath.equities.Calibration.LNSVQDModelCalibrator;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulator;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorMC;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import org.junit.Test;

public class LNSVQDModelCalibratorTest extends TestsSetupForLNSVQD{

	@Test
	public void precalibrateMeanReversion() throws Exception {
		setDAXHestonFebruarySetupSIM();

		/**
		 * 1. Calibrate and get cvalibrated paramerters
		 */
		double[] calibratedParameters;
		double maturity = maturityGrid[maturityGrid.length - 1];
		int numberOfPaths = 10000;
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		LNSVQDPathSimulator lnsvqdPathSimulator = new LNSVQDPathSimulatorMC(valuationDate, disountCurve, equityForwardStructure,
				numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
		calibratedParameters = LNSVQDModelCalibrator.calibrateMRParams(selectedParamsToCalibrate, lnsvqdPathSimulator, acfAtLags);
		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null);
		impliedVolSurface.printVolSurfaceForOutput();
	}
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