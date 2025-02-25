package net.finmath.equities;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.Black76Model;
import net.finmath.equities.models.LNSVQD.*;
import net.finmath.equities.models.VolatilityPointsSurface;
import org.junit.Test;

import java.time.LocalDate;

public class LNSVQDModelCalibratorTest extends TestsSetupForLNSVQD{
	@Test
	public void calibrateTest() throws Exception {
		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(spot0, selectedParams[0], selectedParams[1], selectedParams[2], selectedParams[3], selectedParams[4], selectedParams[5], 0, valuationDate, equityForwardStructure);

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
		calibratedParameters = LNSVQDModelCalibrator.calibrate(selectedParams, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);

		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null);
		impliedVolSurface.printVolSurfaceForOutput();
	}

	@Test
	public void calibrateTestHeston() throws Exception {
		/**
		 * 1. Calibrate and get cvalibrated paramerters
		 */
		double[] calibratedParameters;
		int[] indicesCalibratedParams = {
				/*0
				,*/ 1
				, 2
				/*, 3
				, 4
				, 5*/
		};
		calibratedParameters = LNSVQDModelCalibrator.calibrate(selectedParams, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);

		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null);
		impliedVolSurface.printVolSurfaceForOutput();
	}

	/**
	 * Print calibrated surface
	 */
	@Test
	public void printCalibratedSurfaceMC() throws Exception {
		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(spot0, selectedParams[0], selectedParams[1], selectedParams[2], selectedParams[3], selectedParams[4], selectedParams[5], 0, valuationDate, equityForwardStructure);

		for(VolatilityPoint volPoint : volatilityPointsSurface.getVolatilityPoints()) {
			LocalDate maturity = volPoint.getDate();
			double ttm = dayCountConvention.getDaycountFraction(valuationDate, maturity);
			double strike = volPoint.getStrike();
			double forward = equityForwardStructure.getForward(ttm) * spot0;
			double discFac = equityForwardStructure.getRepoCurve().getDiscountFactor(ttm);
			double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
					ttm, (int) Math.round(ttm * 365.));
			LNSVQDEuropeanPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDEuropeanPriceSimulator(lnsvqdModelAnalyticalPricer, 100000, timeGrid, false);
			lnsvqdCallPriceSimulator.precalculatePaths(35425);
			double price = lnsvqdCallPriceSimulator.getEuropeanPrice(strike, ttm, 1); // TODO: Change!
			double impliedVol = Black76Model.optionImpliedVolatility(forward, strike, ttm, price / discFac, true);
			System.out.println("MC implied vols" );{
				System.out.print(maturity + "\t" + strike + "\t" + impliedVol);
				System.out.print("\n");
			}
		}
	}

}