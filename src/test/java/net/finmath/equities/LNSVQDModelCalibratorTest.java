package net.finmath.equities;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.DynamicVolatilitySurface;
import net.finmath.equities.models.LNSVQD.LNSVQDModel;
import net.finmath.equities.models.LNSVQD.LNSVQDModelAnalyticalPricer;
import net.finmath.equities.models.LNSVQD.LNSVQDModelCalibrator;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.equities.models.VolatilitySurface;
import net.finmath.optimizer.SolverException;
import org.junit.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

public class LNSVQDModelCalibratorTest {
	/**
	 * Model params
	 */
	private final double spot0 = 1;
	private final double sigma0 = 0.5;
	// Value as in paper
	private final double kappa1 = 2.21;
	// Value as in paper
	private final double kappa2 = 2.18;
	private final double theta = 0.5;
	private final double beta = 0.3;
	private final double epsilon = 0.3;

	/**
	 * Get pricer
	 */
	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0);

	/**
	 * Create a strikes * ttm - grid
	 */
	double[] strikes = LNSVQDUtils.createTimeGrid(0.6 * spot0, 1.4 * spot0, 8);
	double[] ttms = LNSVQDUtils.createTimeGrid(1 / 12., 2, 23);

	/**
	 * Other
	 */
	LocalDate today = LocalDate.parse("2024-10-18");

	/**
	 * ***************************************************+
	 * SECTION 1: Get model-implied vol sufrace
	 * ***************************************************+
	 */
	@Test
	public void printModelImpliedVolSurface() {
		VolatilityPoint[] volatilityPoints = new VolatilityPoint[strikes.length * ttms.length];
		for(int i = 0; i < ttms.length; i++) {
			for(int j = 0; j < strikes.length; j++) {
				// LocalDate date = 0; // TODO
				double strike = strikes[j];
				double price = lnsvqdModelAnalyticalPricer.getCallPrice
						(strike, ttms[i], Math.exp(-lnsvqdModelAnalyticalPricer.getRiskFreeRate() * ttms[i]), 0); // TODO: Check
				// volatilityPoints[i * strikes.length + j] = blackScholesOptionImpliedVolatility(...);
			}

		}
		DynamicVolatilitySurface dynamicVolatilitySurface = new DynamicVolatilitySurface(volatilityPoints, today);

		dynamicVolatilitySurface.printVolSurfaceForOutput();
	}



	/**
	 * ***************************************************+
	 * SECTION 2: Calibration test
	 * ***************************************************+
	 */
	@Test
	public void calibrateTest() throws SolverException {
		/**
		 * 1. Create volatility surface
		 * TODO: Create an interface for excel I/O-functionality
		 */
		// Initialize volatilityPoints
		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();

		// Create and adf volatility points
		/*volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.60, 0.6097, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.80, 0.3665, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.00, 0.1374, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.20, 0.212, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.40, 0.319, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.60, 0.476, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.80, 0.3061, lnsvqdModelAnalyticalPricer.getSpot0()));*/
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.00, 0.1508, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.20, 0.1568, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.40, 0.2238, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 0.60, 0.4171, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 0.80, 0.2762, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.00, 0.1493, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.20, 0.1359, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.40, 0.1878, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 0.60, 0.3471, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 0.80, 0.2427, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.00, 0.1511, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.20, 0.1162, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.40, 0.1464, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 0.60, 0.3159, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 0.80, 0.2289, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.00, 0.1545, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.20, 0.1164, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.40, 0.1341, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 0.60, 0.3028, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 0.80, 0.2239, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 1.00, 0.1575, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 1.20, 0.1198, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 1.40, 0.1307, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 0.60, 0.2908, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 0.80, 0.2204, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 1.00, 0.1625, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 1.20, 0.1279, lnsvqdModelAnalyticalPricer.getSpot0()));
		/*volatilityPoints.add(makeVolatilityPoint("2025-12-19", 1.40, 0.1322, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 0.60, 0.2725, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 0.80, 0.2127, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 1.00, 0.1645, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 1.20, 0.1331, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 1.40, 0.1296, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 0.60, 0.2649, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 0.80, 0.2104, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 1.00, 0.167, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 1.20, 0.1383, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 1.40, 0.1321, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 0.60, 0.2584, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 0.80, 0.2087, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 1.00, 0.1694, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 1.20, 0.1424, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 1.40, 0.1333, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 0.60, 0.2544, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 0.80, 0.2079, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.00, 0.1714, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.20, 0.1458, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.40, 0.1351, lnsvqdModelAnalyticalPricer.getSpot0()));*/

		// Calibration day
		LocalDate today = LocalDate.parse("2024-10-18");

		// Create volatility surface
		DynamicVolatilitySurface volatilitySurface = new DynamicVolatilitySurface(volatilityPoints, today);

		/**
		 * 2. ...
		 */

		final double[] initialVolatilityParameters = {kappa1, kappa2, theta, beta, epsilon};

		/**
		 * 3. Calibrate and get cvalibrated paramerters
		 */

		double[] calibratedParameters;
		int[] indicesCalibratedParams = {2, 3, 4};
		calibratedParameters = LNSVQDModelCalibrator.calibrate(initialVolatilityParameters, lnsvqdModelAnalyticalPricer, volatilitySurface, indicesCalibratedParams);


		System.out.println("Calibrated parameters:");
		for(double param : calibratedParameters) {
			System.out.println(param);
		}

	}

	private VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * spot;
		return new VolatilityPoint(maturity, strike, volatility);
	}

}