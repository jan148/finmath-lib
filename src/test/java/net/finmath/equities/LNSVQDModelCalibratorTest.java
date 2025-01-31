package net.finmath.equities;

import net.finmath.equities.marketdata.AffineDividend;
import net.finmath.equities.marketdata.AffineDividendStream;
import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.Black76Model;
import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQD.*;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.optimizer.SolverException;
import net.finmath.time.daycount.DayCountConvention;
import net.finmath.time.daycount.DayCountConvention_ACT_365;
import org.junit.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class LNSVQDModelCalibratorTest {
	/**
	 * Time params
	 */
	LocalDate valuationDate = LocalDate.parse("2024-09-30");
	double spot0 = 19324.93;
	double[] paramVector = new double[]{
			0.16778079819183128,
			10.84439280484746,
			11.081415266404194,
			0.14585436051519188,
			-2.714462050910902,
			2.6723715330684525
	};

	// 0.17971557337087685	2.166285310451572	0.0	0.11259391112260439	-1.1521366420370178	1.211288
	double[] paramVectorHeston = new double[]{
			0.156642865,
			4,
			0,
			0.191510584,
			-0.81402,
			1.211288
	};

	double[] paramVectorHestonNotUsed = new double[]{
			.19095902555549293,
			1.595312176431141,
			0,
			.10219241040931401,
			-1.009187596625914,
			0.9910223497930533
	};

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

		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(spot0, paramVector[0], paramVector[1], paramVector[2], paramVector[3], paramVector[4], paramVector[5], 0, valuationDate, equityForwardStructure);

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
		calibratedParameters = LNSVQDModelCalibrator.calibrate(paramVector, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);

		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilityPointsSurface);
		impliedVolSurface.printVolSurfaceForOutput();
	}

	@Test
	public void calibrateTestHeston() throws Exception {
		setTargetSurfaceHeston();

		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(
						spot0
						, paramVectorHeston[0]
						, paramVectorHeston[1]
						, paramVectorHeston[2]
						, paramVectorHeston[3]
						, paramVectorHeston[4]
						, paramVectorHeston[5]
						, 0, valuationDate, equityForwardStructure);

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
		calibratedParameters = LNSVQDModelCalibrator.calibrate(paramVectorHeston, indicesCalibratedParams, lnsvqdModelAnalyticalPricer, volatilityPointsSurface);

		lnsvqdModelAnalyticalPricer.setVolatilityParameters(calibratedParameters);
		VolatilityPointsSurface impliedVolSurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilityPointsSurface);
		impliedVolSurface.printVolSurfaceForOutput();
	}

	/**
	 * Print calibrated surface
	 */
	@Test
	public void printCalibratedSurface() throws Exception {
		setTargetSurface();

		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(spot0, paramVector[0], paramVector[1], paramVector[2], paramVector[3], paramVector[4], paramVector[5], 0, valuationDate, equityForwardStructure);

		VolatilityPointsSurface modelImpliedVolatilitySurface = lnsvqdModelAnalyticalPricer.getImpliedVolSurface(volatilityPointsSurface);
		modelImpliedVolatilitySurface.printVolSurfaceForOutput();
	}

	@Test
	public void printCalibratedSurfaceMC() throws Exception {
		setTargetSurfaceHeston();

		LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
				new LNSVQDModelAnalyticalPricer(spot0, paramVectorHeston[0], paramVectorHeston[1], paramVectorHeston[2], paramVectorHeston[3], paramVectorHeston[4], paramVectorHeston[5], 0, valuationDate, equityForwardStructure);

		for(VolatilityPoint volPoint : volatilityPointsSurface.getVolatilityPoints()) {
			LocalDate maturity = volPoint.getDate();
			double ttm = dayCountConvention.getDaycountFraction(valuationDate, maturity);
			double strike = volPoint.getStrike();
			double forward = equityForwardStructure.getForward(ttm) * spot0;
			double discFac = equityForwardStructure.getRepoCurve().getDiscountFactor(ttm);
			double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
					ttm, (int) Math.round(ttm * 365.));
			LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, 100000, timeGrid);
			lnsvqdCallPriceSimulator.precalculatePaths(35425);
			double price = lnsvqdCallPriceSimulator.getCallPrice(strike);
			double impliedVol = Black76Model.optionImpliedVolatility(forward, strike, ttm, price / discFac, true);
			System.out.println("MC implied vols" );{
				System.out.print(maturity + "\t" + strike + "\t" + impliedVol);
				System.out.print("\n");
			}
		}
	}

	/**
	 * UTILS
	 */
	private void setTargetSurface() {
		// Initialize volatilityPoints
		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();

		// Create and adf volatility points
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.60, 0.6097, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 0.80, 0.3665, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.00, 0.1374, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.20, 0.212, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-10-18", 1.40, 0.319, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.60, 0.476, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 0.80, 0.3061, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.00, 0.1508, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.20, 0.1568, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-11-15", 1.40, 0.2238, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 0.60, 0.4171, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 0.80, 0.2762, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.00, 0.1493, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.20, 0.1359, spot0));
		volatilityPoints.add(makeVolatilityPoint("2024-12-20", 1.40, 0.1878, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 0.60, 0.3471, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 0.80, 0.2427, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.00, 0.1511, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.20, 0.1162, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-03-21", 1.40, 0.1464, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 0.60, 0.3159, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 0.80, 0.2289, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.00, 0.1545, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.20, 0.1164, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-06-20", 1.40, 0.1341, spot0));
		/*volatilityPoints.add(makeVolatilityPoint("2025-09-19", 0.60, 0.3028, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 0.80, 0.2239, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 1.00, 0.1575, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 1.20, 0.1198, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-09-19", 1.40, 0.1307, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 0.60, 0.2908, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 0.80, 0.2204, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 1.00, 0.1625, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 1.20, 0.1279, spot0));
		volatilityPoints.add(makeVolatilityPoint("2025-12-19", 1.40, 0.1322, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 0.60, 0.2725, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 0.80, 0.2127, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 1.00, 0.1645, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 1.20, 0.1331, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-06-19", 1.40, 0.1296, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 0.60, 0.2649, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 0.80, 0.2104, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 1.00, 0.167, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 1.20, 0.1383, spot0));
		volatilityPoints.add(makeVolatilityPoint("2026-12-18", 1.40, 0.1321, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 0.60, 0.2584, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 0.80, 0.2087, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 1.00, 0.1694, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 1.20, 0.1424, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-06-18", 1.40, 0.1333, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 0.60, 0.2544, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 0.80, 0.2079, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.00, 0.1714, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.20, 0.1458, spot0));
		volatilityPoints.add(makeVolatilityPoint("2027-12-17", 1.40, 0.1351, spot0));*/

		// Create volatility surface
		volatilityPointsSurface = new VolatilityPointsSurface(volatilityPoints, valuationDate, dayCountConvention);
	}


	private void setTargetSurfaceHeston() {
		double[] ttms = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		LocalDate[] dates = Arrays.stream(ttms).mapToObj(ttm -> {
					long days = Math.round(ttm * 365);
					return valuationDate.plusDays(days);}).toArray(LocalDate[]::new);

		// Initialize volatilityPoints
		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();
		/*0.401531086	0.262002147	0.147047443	0.131142026	0.182123154
		0.33754687	0.23199757	0.145062873	0.117617253	0.147591954
		0.30729848	0.218394871	0.146522784	0.117841898	0.137581478
		0.293414523	0.213139202	0.149207569	0.120862294	0.134960126
		0.280989633	0.20937458	0.15356843	0.126817863	0.135592582
		0.269761611	0.204298557	0.153865947	0.128339944	0.133914086*/

		// Create and adf volatility points
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.60, 0.401531086, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.80, 0.262002147, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.00, 0.147047443, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.20, 0.131142026, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.40, 0.182123154, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.60, 0.33754687, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.80, 0.23199757, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.00, 0.145062873, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.20, 0.117617253, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.40, 0.147591954, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.60, 0.30729848, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.80, 0.218394871, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.00, 0.146522784, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.20, 0.117841898, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.40, 0.137581478, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.60, 0.293414523, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.80, 0.213139202, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.00, 0.149207569, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.20, 0.120862294, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.40, 0.134960126, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.60, 0.280989633, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.80, 0.20937458, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.00, 0.15356843, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.20, 0.126817863, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.40, 0.135592582, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.60, 0.269761611, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.80, 0.204298557, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.00, 0.153865947, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.20, 0.128339944, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.40, 0.133914086, spot0));


		// Create volatility surface
		volatilityPointsSurface = new VolatilityPointsSurface(volatilityPoints, valuationDate, dayCountConvention);
	}

	private VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * spot * equityForwardStructure.getForward(maturity);
		return new VolatilityPoint(maturity, strike, volatility);
	}

}