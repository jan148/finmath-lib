package net.finmath.equities;

import net.finmath.equities.marketdata.*;
import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQD.LNSVQDEuropeanPriceSimulator;
import net.finmath.equities.models.LNSVQD.LNSVQDModelAnalyticalPricer;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.time.daycount.DayCountConvention;
import net.finmath.time.daycount.DayCountConvention_ACT_365;
import org.apache.commons.math3.util.Pair;
import org.junit.jupiter.api.BeforeEach;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public abstract class TestsSetupForLNSVQD {
	/**
	 * Time params
	 */
	LocalDate valuationDate = LocalDate.parse("2024-09-30");
	double spot0 = 1;

	/**
	 * Create analytical model
	 */
	final double[] paramVectorBitcoin = new double[]{
			0.8376,
			3.1844,
			3.058,
			1.0413,
			0.1514,
			1.8458
	};
	final double[] paramBlackScholes = new double[]{
			0.2,
			0,
			0,
			0,
			0,
			0
	};

	final double[] paramVectorLowVolOfVol = new double[]{
			0.156642865,
			4,
			0,
			0.191510584,
			-0.2,
			0.211288
	};

	final double[] paramVectorInitial = new double[]{
			0.156642865,
			4,
			0,
			0.191510584,
			-0.81402,
			1.211288
	};

	final double[] paramVectorCalibrated = new double[]{
			0.13731067,
			4,
			0,
			0.14664391,
			-1.4890642,
			1.53513925
	};

	final double[] paramVectorCalibratedWPositiveKappa2 = new double[]{
			0.13731067,
			4,
			0.4,
			0.14664391,
			-1.4890642,
			1.53513925
	};

	double[] selectedParams = null;

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
	EquityForwardStructure equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, spot0, yieldCurve, affineDividendStream, dayCountConvention);

	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer = null;

	/**
	 * Create simulation model (not finmath)
	 */
	int numberOfPaths = 100000;
	// TODO: Discounts dates should be decoupled from maturities
	/*double[] maturityGrid = Arrays.stream(discountDates)
			.mapToDouble(date -> dayCountConvention.getDaycountFraction(valuationDate, date))
			.toArray();*/

	/**
	 * Declare volatility surface; will be insatntiated later
	 */
	double[] maturityGrid = null;
	VolatilityPointsSurface volatilityPointsSurface = null;

	/**
	 * ***************************************************+
	 * UTILS
	 * ***************************************************+
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


	public ArrayList<Pair<Double, Double>> setDAXHestonSetup() {
		selectedParams = paramVectorCalibrated;

		maturityGrid = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		double[] forwards = Arrays.stream(maturityGrid).map(x -> equityForwardStructure.getForward(x)).toArray();

		lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(
				spot0
				, selectedParams[0]
				, selectedParams[1]
				, selectedParams[2]
				, selectedParams[3]
				, selectedParams[4]
				, selectedParams[5]
				, 0, valuationDate, equityForwardStructure);

		// Initialize volatilityPoints
		ArrayList<Pair<Double, Double>> strikeMatPairs = new ArrayList<>();

		// Create and adf volatility points
		strikeMatPairs.add(new Pair<>(maturityGrid[0], 0.60 * forwards[0]));
		strikeMatPairs.add(new Pair<>(maturityGrid[0], 0.80 * forwards[0]));
		strikeMatPairs.add(new Pair<>(maturityGrid[0], 1.00 * forwards[0]));
		strikeMatPairs.add(new Pair<>(maturityGrid[0], 1.20 * forwards[0]));
		strikeMatPairs.add(new Pair<>(maturityGrid[0], 1.40 * forwards[0]));
		strikeMatPairs.add(new Pair<>(maturityGrid[1], 0.60 * forwards[1]));
		strikeMatPairs.add(new Pair<>(maturityGrid[1], 0.80 * forwards[1]));
		strikeMatPairs.add(new Pair<>(maturityGrid[1], 1.00 * forwards[1]));
		strikeMatPairs.add(new Pair<>(maturityGrid[1], 1.20 * forwards[1]));
		strikeMatPairs.add(new Pair<>(maturityGrid[1], 1.40 * forwards[1]));
		strikeMatPairs.add(new Pair<>(maturityGrid[2], 0.60 * forwards[2]));
		strikeMatPairs.add(new Pair<>(maturityGrid[2], 0.80 * forwards[2]));
		strikeMatPairs.add(new Pair<>(maturityGrid[2], 1.00 * forwards[2]));
		strikeMatPairs.add(new Pair<>(maturityGrid[2], 1.20 * forwards[2]));
		strikeMatPairs.add(new Pair<>(maturityGrid[2], 1.40 * forwards[2]));
		strikeMatPairs.add(new Pair<>(maturityGrid[3], 0.60 * forwards[3]));
		strikeMatPairs.add(new Pair<>(maturityGrid[3], 0.80 * forwards[3]));
		strikeMatPairs.add(new Pair<>(maturityGrid[3], 1.00 * forwards[3]));
		strikeMatPairs.add(new Pair<>(maturityGrid[3], 1.20 * forwards[3]));
		strikeMatPairs.add(new Pair<>(maturityGrid[3], 1.40 * forwards[3]));
		strikeMatPairs.add(new Pair<>(maturityGrid[4], 0.60 * forwards[4]));
		strikeMatPairs.add(new Pair<>(maturityGrid[4], 0.80 * forwards[4]));
		strikeMatPairs.add(new Pair<>(maturityGrid[4], 1.00 * forwards[4]));
		strikeMatPairs.add(new Pair<>(maturityGrid[4], 1.20 * forwards[4]));
		strikeMatPairs.add(new Pair<>(maturityGrid[4], 1.40 * forwards[4]));
		strikeMatPairs.add(new Pair<>(maturityGrid[5], 0.60 * forwards[5]));
		strikeMatPairs.add(new Pair<>(maturityGrid[5], 0.80 * forwards[5]));
		strikeMatPairs.add(new Pair<>(maturityGrid[5], 1.00 * forwards[5]));
		strikeMatPairs.add(new Pair<>(maturityGrid[5], 1.20 * forwards[5]));
		strikeMatPairs.add(new Pair<>(maturityGrid[5], 1.40 * forwards[5]));

		return strikeMatPairs;
	}

	public ArrayList<Pair<Double, Double>> setBTCSetup() {
		double spot0 = 67843.219;
		selectedParams = paramVectorBitcoin;

		yieldCurve = new FlatYieldCurve(valuationDate, 0., dayCountConvention);

		affineDividends = new AffineDividend[]{new AffineDividend(valuationDate, 0., 0.)};
		affineDividendStream = new AffineDividendStream(affineDividends);
		equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, spot0, yieldCurve, affineDividendStream, dayCountConvention);

		lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(
						spot0
						, selectedParams[0]
						, selectedParams[1]
						, selectedParams[2]
						, selectedParams[3]
						, selectedParams[4]
						, selectedParams[5]
						, 0, valuationDate, equityForwardStructure);

		double ttm = 0.10122575874485597;
		maturityGrid = new double[]{0.10122575874485597};

		// Initialize volatilityPoints
		ArrayList<Pair<Double, Double>> strikeMatPairs = new ArrayList<>();

		// Create and adf volatility points
		strikeMatPairs.add(new Pair(ttm, 45000.));
		strikeMatPairs.add(new Pair(ttm, 48000.));
		strikeMatPairs.add(new Pair(ttm, 55000.));
		strikeMatPairs.add(new Pair(ttm, 58000.));
		strikeMatPairs.add(new Pair(ttm, 64000.));
		strikeMatPairs.add(new Pair(ttm, 65000.));
		strikeMatPairs.add(new Pair(ttm, 70000.));
		strikeMatPairs.add(new Pair(ttm, 75000.));
		strikeMatPairs.add(new Pair(ttm, 80000.));
		strikeMatPairs.add(new Pair(ttm, 85000.));
		strikeMatPairs.add(new Pair(ttm, 90000.));
		strikeMatPairs.add(new Pair(ttm, 10000.));
		strikeMatPairs.add(new Pair(ttm, 12000.));

		return strikeMatPairs;
	}

	private VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * spot * equityForwardStructure.getForward(maturity);
		return new VolatilityPoint(maturity, strike, volatility);
	}

}
