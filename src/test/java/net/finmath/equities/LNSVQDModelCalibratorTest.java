package net.finmath.equities;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.equities.models.LNSVQD.LNSVQDModelAnalyticalPricer;
import net.finmath.equities.models.LNSVQD.LNSVQDModelCalibrator;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.optimizer.SolverException;
import net.finmath.time.daycount.DayCountConvention;
import net.finmath.time.daycount.DayCountConvention_ACT_365;
import org.junit.Test;

import java.time.LocalDate;
import java.util.ArrayList;

public class LNSVQDModelCalibratorTest {
	/**
	 * Time params
	 */
	LocalDate valuationDate = LocalDate.parse("2024-09-18");

	/**
	 * Model params
	 */
	// Right params: sigma0=0.8327, theta=1.0139, kappa1=4.8606, kappa2=4.7938, beta=0.1985, volvol=2.3690
	// sigma0=0.8376, theta=1.0413, kappa1=3.1844, kappa2=3.058, beta=0.1514, volvol=1.8458
	private final double spot0 = 1;
	private final double sigma0 = 0.1;
	// Value as in paper
	private final double kappa1 = 0;
	// Value as in paper
	private final double kappa2 = 0;
	private final double theta = 0;
	private final double beta = 0;
	private final double epsilon = 0;

	/**
	 * Get pricer
	 */
	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0, valuationDate);

	/**
	 * Create a strikes * ttm - grid
	 */
	double[] strikes = LNSVQDUtils.createTimeGrid(0.6 * spot0, 1.4 * spot0, 4);
	double[] ttms = LNSVQDUtils.createTimeGrid(1 / 12., 2, 23);

	/**
	 * Create volatilty surface
	 */
	ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();
	VolatilityPointsSurface volatilityPointsSurface;

	/**
	 * Other
	 */
	DayCountConvention dayCountConvention = new DayCountConvention_ACT_365();
	int daysPerYear = 365;

	/**
	 * Other
	 */

	/**
	 * ***************************************************+
	 * SECTION 1: Test model-implied vol-surface
	 * ***************************************************+
	 */
	@Test
	public void getImpliedVolSurfaceTest() {

	}

	/**
	 * ***************************************************+
	 * SECTION 2: Get model-implied vol sufrace
	 * ***************************************************+
	 */
	@Test
	public void printModelImpliedVolSurface() {
		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();
		for(int i = 0; i < ttms.length; i++) {
			double ttm = ttms[i];
			double discountFactor = Math.exp(-lnsvqdModelAnalyticalPricer.getRiskFreeRate(ttm) * ttm);
			double forward = lnsvqdModelAnalyticalPricer.getSpot0() / discountFactor;
			for(int j = 0; j < strikes.length; j++) {
				LocalDate date = valuationDate.plusDays(Math.round(daysPerYear * ttm));
				double daysInBetween = dayCountConvention.getDaycount(valuationDate, date);
				double strike = strikes[j];
				double price = lnsvqdModelAnalyticalPricer.getCallPrice
						(strike, ttms[i], discountFactor, 0); // TODO: Check
				/*double price = AnalyticFormulas.blackScholesOptionValue(spot0, lnsvqdModelAnalyticalPricer.getRiskFreeRate(), sigma0, ttm, strike);*/
				double impliedVol = AnalyticFormulas.blackScholesOptionImpliedVolatility
						(forward, ttm, strike, discountFactor, price);
				volatilityPoints.add(new VolatilityPoint(date, strike, impliedVol));
				// System.out.println(impliedVol);
			}
		}
		VolatilityPointsSurface volatilityPointsSurface = new VolatilityPointsSurface(volatilityPoints, valuationDate, dayCountConvention);
		volatilityPointsSurface.printVolSurfaceForOutput();
	}

	/**
	 * ***************************************************+
	 * SECTION 2: Calibration test
	 * ***************************************************+
	 */
	@Test
	public void calibrateTest() throws SolverException {
		setTargetSurface();

		/**
		 * 0. ...
		 */
		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(spot0, sigma0, kappa1, kappa2, theta, beta, epsilon, 0, valuationDate);

		/**
		 * 1. ...
		 */
		final double[] initialVolatilityParameters = {sigma0, kappa1, kappa2, theta, beta, epsilon};

		/**
		 * 2. Calibrate and get cvalibrated paramerters
		 */
		double[] calibratedParameters;
		int[] indicesCalibratedParams = {0, 1, 2 , 3, 4, 5};
		calibratedParameters = LNSVQDModelCalibrator.calibrate(initialVolatilityParameters, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);

		System.out.println("Calibrated parameters:");
		for(double param : calibratedParameters) {
			System.out.println(param);
		}

		try {
			VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilityPointsSurface);
			impliedVolSurface.printVolSurfaceForOutput();
		} catch(Exception e) {
			throw new RuntimeException(e);
		}


	}

	/**
	 * Print calibrated surface
	 */
	@Test
	public void printCalibratedSurface() throws Exception {
		setTargetSurface();
		VolatilityPointsSurface modelImpliedVolatilitySurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilityPointsSurface);
		modelImpliedVolatilitySurface.printVolSurfaceForOutput();
	}

	/**
	 * UTILS
	 */
	private void setTargetSurface() {
		// Initialize volatilityPoints
		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();

		// Create and adf volatility points
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.60, 0.6097, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.80, 0.3665, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.00, 0.1374, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.20, 0.212, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.40, 0.319, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.60, 0.476, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.80, 0.3061, lnsvqdModelAnalyticalPricer.getSpot0()));
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
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 1.40, 0.1322, lnsvqdModelAnalyticalPricer.getSpot0()));
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
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.40, 0.1351, lnsvqdModelAnalyticalPricer.getSpot0()));

		// Create volatility surface
		volatilityPointsSurface = new VolatilityPointsSurface(volatilityPoints, valuationDate, dayCountConvention);
	}

	private VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * spot;
		return new VolatilityPoint(maturity, strike, volatility);
	}

}