package net.finmath.equities;

import net.finmath.equities.marketdata.AffineDividend;
import net.finmath.equities.marketdata.AffineDividendStream;
import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQD.LNSVQDModel;
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
import java.util.Random;

public class LNSVQDModelCalibratorTest {
	/**
	 * Time params
	 */
	LocalDate valuationDate = LocalDate.parse("2024-09-30");
	double spot0 = 1;

	/**
	 * Create volatilty surface
	 */
	VolatilityPointsSurface volatilityPointsSurface;

	/**
	 * Other
	 */
	DayCountConvention dayCountConvention = new DayCountConvention_ACT_365();

	/**
	 * Create forward stucture
	 */
	LocalDate[] discountDates = LNSVQDUtils.createLocalDateList(new String[]{
			"2024-10-18",
			"2024-11-15",
			"2024-12-20",
			"2025-03-21",
			"2025-06-20",
			"2025-09-19",
			"2025-12-19",
			"2026-06-19",
			"2026-12-18",
			"2027-06-18",
			"2027-12-17"
	});
	double[] discountFactors = new double[]{
			0.997589468
			, 0.994816096
			, 0.991454403
			, 0.982961371
			, 0.973114616
			, 0.96776706
			, 0.962814667
			, 0.949011402
			, 0.938160121
			, 0.923490121
			, 0.91252815
	};
	YieldCurve yieldCurve = new YieldCurve("Discount curve"
			, valuationDate
			, dayCountConvention
			, discountDates
			, discountFactors);

	AffineDividend[] affineDividends = new AffineDividend[]{new AffineDividend(valuationDate, 0., 0.)};
	AffineDividendStream affineDividendStream = new AffineDividendStream(affineDividends);
	EquityForwardStructure equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, 1, yieldCurve, affineDividendStream, dayCountConvention);

	/**
	 * ***************************************************+
	 * SECTION 1: Calibration test
	 * ***************************************************+
	 */
	@Test
	public void calibrateTest() throws Exception {
		setTargetSurface();
		Random r = new Random();
		double lowestError = 100000000;
		double[] bestParam;

		double[] paramVector = new double[]{
				0.20216505228211554
				, 2.7660217303300856
				, 16.639965808971205
				, 0.2360647350930049
				, -1.693221137952367
				, 0.9430592760276209};

		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(spot0, paramVector[0], paramVector[1], paramVector[2], paramVector[3], paramVector[4], paramVector[5], 0, valuationDate, equityForwardStructure);

		/**
		 * 1. Calibrate and get cvalibrated paramerters
		 */
		double[] calibratedParameters;
		int[] indicesCalibratedParams = {0, /*1, 2,*/ 3, 4, 5};
		calibratedParameters = LNSVQDModelCalibrator.calibrate(paramVector, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);

		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilityPointsSurface);
		impliedVolSurface.printVolSurfaceForOutput();
	}

	/**
	 * Print calibrated surface
	 */
	/*@Test
	public void printCalibratedSurface() throws Exception {
		setTargetSurface();
		VolatilityPointsSurface modelImpliedVolatilitySurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilityPointsSurface);
		modelImpliedVolatilitySurface.printVolSurfaceForOutput();
	}*/

	/**
	 * UTILS
	 */
	private void setTargetSurface() {
		// Initialize volatilityPoints
		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();

		// Create and adf volatility points
		/*volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.60, 0.6097, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.80, 0.3665, lnsvqdModelAnalyticalPricer.getSpot0()));*/
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.00, 0.1374, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.20, 0.212, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.40, 0.319, spot0));
		/*volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.60, 0.476, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.80, 0.3061, lnsvqdModelAnalyticalPricer.getSpot0()));*/
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.00, 0.1508, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.20, 0.1568, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.40, 0.2238, spot0));
		/*volatilityPoints.add(makeVolatilityPoint("2024-12-20", 0.60, 0.4171, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 0.80, 0.2762, lnsvqdModelAnalyticalPricer.getSpot0()));*/
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.00, 0.1493, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.20, 0.1359, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.40, 0.1878, spot0));
		/*volatilityPoints.add(makeVolatilityPoint("2025-03-21", 0.60, 0.3471, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 0.80, 0.2427, lnsvqdModelAnalyticalPricer.getSpot0()));*/
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.00, 0.1511, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.20, 0.1162, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.40, 0.1464, spot0));
		/*volatilityPoints.add(makeVolatilityPoint("2025-06-20", 0.60, 0.3159, lnsvqdModelAnalyticalPricer.getSpot0()));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 0.80, 0.2289, lnsvqdModelAnalyticalPricer.getSpot0()));*/
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.00, 0.1545, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.20, 0.1164, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.40, 0.1341, spot0));
		/*volatilityPoints.add(makeVolatilityPoint("2025-09-19", 0.60, 0.3028, lnsvqdModelAnalyticalPricer.getSpot0()));
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
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.40, 0.1351, lnsvqdModelAnalyticalPricer.getSpot0()));*/

		// Create volatility surface
		volatilityPointsSurface = new VolatilityPointsSurface(volatilityPoints, valuationDate, dayCountConvention);
	}

	private VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * spot;
		return new VolatilityPoint(maturity, strike, volatility);
	}

}