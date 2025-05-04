package net.finmath.equities;

import net.finmath.equities.Calibration.LNSVQDModelCalibrator;
import net.finmath.equities.Simulation.PathSimulator;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.equities.pricer.LNSVQDModelAnalyticalPricer;
import org.junit.Test;

public class LNSVQDModelCalibratorTest extends TestsSetupForLNSVQD{

	@Test
	public void printModelACF() throws Exception {
		loadS20(); // loadS25();

		int startingIndex = 1;
		double[] startingValueLNSVQD = new double[]{spot0, selectedParamsLNSVQD[0]};
		double maturity = 10; //maturityGrid[maturityGrid.length - 1];
		int numberOfPaths = 100000;
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(new double[]{},
						(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		PathSimulator lnsvqdPathSimulator = new PathSimulator(valuationDate, disountCurve, equityForwardStructure,
				numberOfPaths, timeGrid, maturityGrid);
		lnsvqdPathSimulator.precalculatePaths(3105, false, startingIndex, startingValueLNSVQD, Boolean.TRUE
				, "LnsvqdForwardEuler", "MC", lnsvqdModelDescriptor, null);
		double[][] volaPaths = lnsvqdPathSimulator.path[1];
		double[] acs = LNSVQDModelCalibrator.getAcfsFromVolPaths(volaPaths, acfAtLags.size(), 10);
		LNSVQDUtils.printArray(acs);
	}

	@Test
	public void calibrateTest() throws Exception {
		loadS20();

		double[] calibratedParameters;
		int[] indicesCalibratedParams = {
				0
				/*, 1
				, 2*/
				, 3
				, 4
				, 5
		};
		calibratedParameters = LNSVQDModelCalibrator.calibrate(selectedParamsToCalibrate, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);
		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricerCalib = lnsvqdModelAnalyticalPricer.copyWithNewParameters(
				lnsvqdModelAnalyticalPricer.getSpot0(), calibratedParameters[0], calibratedParameters[1], calibratedParameters[2], calibratedParameters[3]
				, calibratedParameters[4], calibratedParameters[5]
				, lnsvqdModelAnalyticalPricer.getI0(), lnsvqdModelAnalyticalPricer.getSpotDate(), lnsvqdModelAnalyticalPricer.getDiscountCurve(), lnsvqdModelAnalyticalPricer.getEquityForwardStructure());
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricerCalib.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null);
		impliedVolSurface.printVolSurfaceForOutput();
	}

}