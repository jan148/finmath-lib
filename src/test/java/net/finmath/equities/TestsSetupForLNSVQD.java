package net.finmath.equities;

import net.finmath.equities.marketdata.*;
import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.pricer.LNSVQDModelAnalyticalPricer;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.time.daycount.DayCountConvention;
import net.finmath.time.daycount.DayCountConvention_ACT_365;
import org.apache.commons.math3.util.Pair;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;

public abstract class TestsSetupForLNSVQD {
	/**
	 * Time params
	 */
	static LocalDate valuationDate = LocalDate.parse("2024-09-30");
	static double spot0 = 1;

	/**
	 * Create analytical model
	 */
	static double[] paramVectorBitcoin = new double[]{
			0.8376,
			3.1844,
			3.058,
			1.0413,
			0.1514,
			1.8458
	};

	static double[] paramBlackScholes = new double[]{
			0.2,
			0,
			0,
			0,
			0,
			0
	};

	static double[] paramVectorLowVolOfVol = new double[]{
			0.156642865,
			4,
			0,
			0.191510584,
			-0.2,
			0.211288
	};

	static double[] paramVectorInitial = new double[]{
			0.156642865,
			4,
			0,
			0.191510584,
			-0.81402,
			1.211288
	};

	static double[] paramVectorCalibrated = new double[]{
			0.13731016014532768,
			4,
			0,
			0.14664386051912978,
			-1.4890940544243099,
			1.5351312426087254
	};

	static double[] paramVectorCalibratedWPositiveKappa2 = new double[]{
			0.13731067,
			4,
			0.4,
			0.14664391,
			-1.4890642,
			1.53513925
	};

	static double[] paramVectorHeston = new double[]{
			0.024536987
			, 4.
			, 0.036676304
			, 1.211288333
			, -0.672024524
	};

	static double[] selectedParamsLNSVQD;
	static double[] selectedParamsHeston;

	/**
	 * Other
	 */
	static DayCountConvention dayCountConvention = new DayCountConvention_ACT_365();

	/**
	 * Create forward stucture
	 */
	static LocalDate[] discountDates;

	static double[] discountFactors;

	static YieldCurve disountCurve;
	static LocalDate[] forwardDates;

	static double[] forwardDiscountFactors;

	static YieldCurve forwardCurve;

	static AffineDividend[] affineDividends;
	static AffineDividendStream affineDividendStream;
	static EquityForwardStructure equityForwardStructure;

	static LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer;

	/**
	 * Create simulation model (not finmath)
	 */
	static int numberOfPaths = 100000;

	/**
	 * Declare volatility surface; will be insatntiated later
	 */
	static double[] maturityGrid;
	static VolatilityPointsSurface volatilityPointsSurface;

	/**
	 * ***************************************************+
	 * UTILS
	 * ***************************************************+
	 */
	public static void setDAXHestonSetupCALIB() {
		spot0 = 1;
		selectedParamsLNSVQD = paramVectorInitial;
		selectedParamsHeston = paramVectorHeston;

		double[] ttms = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		LocalDate[] dates = Arrays.stream(ttms).mapToObj(ttm -> {
			long days = Math.round(ttm * 365);
			return valuationDate.plusDays(days);
		}).toArray(LocalDate[]::new);

		discountDates = LNSVQDUtils.createLocalDateList(new String[]{
				"2024-10-07",
				"2024-10-14",
				"2024-10-30",
				"2024-11-29",
				"2024-12-30",
				"2025-01-30",
				"2025-02-28",
				"2025-03-31",
				"2025-04-30",
				"2025-05-30",
				"2025-06-30",
				"2025-07-30",
				"2025-08-29",
				"2025-09-30",
				"2026-03-30",
				"2026-09-30",
				"2027-09-30",
				"2028-09-29"
		});

		discountFactors = new double[]{
				0.99933661
				, 0.9986739
				, 0.99722853
				, 0.99457625
				, 0.99197353
				, 0.98951738
				, 0.98740534
				, 0.98533962
				, 0.9834937
				, 0.98180536
				, 0.98010916
				, 0.97865945
				, 0.97716886
				, 0.97558229
				, 0.96739905
				, 0.95894793
				, 0.94120254
				, 0.92255552
		};

		disountCurve = new YieldCurve("Discount curve"
				, valuationDate
				, dayCountConvention
				, discountDates
				, discountFactors);

		forwardDates = LNSVQDUtils.createLocalDateList(new String[]{
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

		forwardDiscountFactors = new double[]{
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

		forwardCurve = new YieldCurve("Discount curve"
				, valuationDate
				, dayCountConvention
				, forwardDates
				, forwardDiscountFactors);

		affineDividends = new AffineDividend[]{new AffineDividend(valuationDate, 0., 0.)};
		affineDividendStream = new AffineDividendStream(affineDividends);
		equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, spot0, forwardCurve, affineDividendStream, dayCountConvention);

		lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(
				spot0
				, selectedParamsLNSVQD[0]
				, selectedParamsLNSVQD[1]
				, selectedParamsLNSVQD[2]
				, selectedParamsLNSVQD[3]
				, selectedParamsLNSVQD[4]
				, selectedParamsLNSVQD[5]
				, 0
				, valuationDate
				, disountCurve
				, equityForwardStructure);

		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<VolatilityPoint>();

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


	public ArrayList<Pair<Double, Double>> setDAXHestonSetupSIM() {
		spot0 = 1;
		selectedParamsLNSVQD = paramVectorCalibrated;
		selectedParamsHeston = paramVectorHeston;

		double[] ttms = new double[]{0.25/*, 0.5, 0.75, 1, 1.25, 1.5*/};
		LocalDate[] dates = Arrays.stream(ttms).mapToObj(ttm -> {
			long days = Math.round(ttm * 365);
			return valuationDate.plusDays(days);
		}).toArray(LocalDate[]::new);

		discountDates = LNSVQDUtils.createLocalDateList(new String[]{
				"2024-10-07",
				"2024-10-14",
				"2024-10-30",
				"2024-11-29",
				"2024-12-30",
				"2025-01-30",
				"2025-02-28",
				"2025-03-31",
				"2025-04-30",
				"2025-05-30",
				"2025-06-30",
				"2025-07-30",
				"2025-08-29",
				"2025-09-30",
				"2026-03-30",
				"2026-09-30",
				"2027-09-30",
				"2028-09-29"
		});

		discountFactors = new double[]{
				0.99933661
				, 0.9986739
				, 0.99722853
				, 0.99457625
				, 0.99197353
				, 0.98951738
				, 0.98740534
				, 0.98533962
				, 0.9834937
				, 0.98180536
				, 0.98010916
				, 0.97865945
				, 0.97716886
				, 0.97558229
				, 0.96739905
				, 0.95894793
				, 0.94120254
				, 0.92255552
		};

		disountCurve = new YieldCurve("Discount curve"
				, valuationDate
				, dayCountConvention
				, discountDates
				, discountFactors);


		forwardDates = LNSVQDUtils.createLocalDateList(new String[]{
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

		forwardDiscountFactors = new double[]{
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

		forwardCurve = new YieldCurve("Discount curve"
				, valuationDate
				, dayCountConvention
				, forwardDates
				, forwardDiscountFactors);

		affineDividends = new AffineDividend[]{new AffineDividend(valuationDate, 0., 0.)};
		affineDividendStream = new AffineDividendStream(affineDividends);
		equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, spot0, forwardCurve, affineDividendStream, dayCountConvention);

		lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(
				spot0
				, selectedParamsLNSVQD[0]
				, selectedParamsLNSVQD[1]
				, selectedParamsLNSVQD[2]
				, selectedParamsLNSVQD[3]
				, selectedParamsLNSVQD[4]
				, selectedParamsLNSVQD[5]
				, 0
				, valuationDate
				, disountCurve
				, equityForwardStructure);

		maturityGrid = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		double[] forwards = Arrays.stream(maturityGrid).map(x -> equityForwardStructure.getForward(x)).toArray();

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

	public ArrayList<Pair<Double, Double>> setBTCSetupSIM() {
		spot0 = 67843.219;
		selectedParamsLNSVQD = paramVectorBitcoin;
		selectedParamsHeston = paramVectorHeston;

		disountCurve = new FlatYieldCurve(valuationDate, 0., dayCountConvention);
		forwardCurve = new FlatYieldCurve(valuationDate, 0., dayCountConvention);

		affineDividends = new AffineDividend[]{new AffineDividend(valuationDate, 0., 0.)};
		affineDividendStream = new AffineDividendStream(affineDividends);
		equityForwardStructure = new BuehlerDividendForwardStructure(valuationDate, spot0, forwardCurve, affineDividendStream, dayCountConvention);

		lnsvqdModelAnalyticalPricer = new LNSVQDModelAnalyticalPricer(
				spot0
				, selectedParamsLNSVQD[0]
				, selectedParamsLNSVQD[1]
				, selectedParamsLNSVQD[2]
				, selectedParamsLNSVQD[3]
				, selectedParamsLNSVQD[4]
				, selectedParamsLNSVQD[5]
				, 0
				, valuationDate
				, disountCurve
				, equityForwardStructure);

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

	static VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * spot * equityForwardStructure.getForward(maturity);
		return new VolatilityPoint(maturity, strike, volatility);
	}

}
