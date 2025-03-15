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
import java.util.List;

public abstract class TestsSetupForLNSVQD {
	/**
	 * Time params
	 */
	static LocalDate valuationDate;
	static double spot0;

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

	static double[] paramVectorInitial = new double[]{
			0.156642865,
			4,
			0,
			0.191510584,
			-0.81402,
			1.211288
	};

	static double[] paramVectorInitialWithPreCalibratedKappas = new double[]{
			0.156642865,
			4.012705001125997,
			4,
			0.19481803777387482,
			0,
			0.8975621744372831
	};

	static double[] paramVectorCalibrated = new double[]{
			0.13731016014532768,
			4,
			0,
			0.14664386051912978,
			-1.4890940544243099,
			1.5351312426087254
	};

	static double[] paramVectorInitialMarch = new double[]{
			0.768464,
			4,
			0,
			0.497749705,
			-0.81402,
			1.211288
	};

	static double[] paramVectorInitialMarchCalibrated = new double[]{
			0.7620015635018382,
			4.0,
			0.0,
			0.47376982758578084,
			-0.6857270328141885,
			0.8741962423663083
	};

	static double[] paramVectorInitialFebruary = new double[]{
			0.1598567708165421,
			4,
			0,
			0.189967208,
			-0.40790656321439606,
			0.20752601856648067
	};

	static double[] paramVectorInitialFebruaryCalibrated = new double[]{
			0.1598580238458368,
			4,
			0,
			0.1436794749651145,
			-1.452186033001682,
			1.6262530348098245
	};

	static double[] paramVectorInitialMarchCalibratedImproved = new double[]{
			0.6921704766354998,
			0.5160497941112346,
			0.0,
			0.19543814508348928,
			-0.40790656321439606,
			0.20752601856648067
	};

	static double[] paramVectorHeston = new double[]{
			0.024536987
			, 4.
			, 0.036676304
			, 1.211288333
			, -0.672024524
	};

	static double[] paramVectorHestonMarch = new double[]{
			0.590537467
			, 4
			, 0.247754769
			, 0.813577246
			, -0.994513338
	};

	static double[] paramVectorHestonFebruary = new double[]{
			0.034583515
			, 4
			, 0.03608754
			, 1.2835035
			, -0.634272101
	};

	static double[] selectedParamsLNSVQD;
	static double[] selectedParamsHeston;
	static double[] selectedParamsToCalibrate;

	/**
	 * Other
	 */
	static DayCountConvention dayCountConvention = new DayCountConvention_ACT_365();
	static List<Pair<Double, Double>> acfAtLags;

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
	static int numberOfPaths = 100;

	/**
	 * Declare volatility surface; will be insatntiated later
	 */
	static double[] maturityGrid;
	static VolatilityPointsSurface volatilityPointsSurface;
	static ArrayList<Pair<Double, Double>> strikeMatPairs;

	/**
	 * ***************************************************+
	 * UTILS
	 * ***************************************************+
	 */
	public void setDAXHestonSetupSIM() {
		valuationDate = LocalDate.parse("2024-09-30");
		spot0 = 1;
		selectedParamsLNSVQD = paramVectorCalibrated;
		selectedParamsHeston = paramVectorHeston;
		selectedParamsToCalibrate = paramVectorInitialWithPreCalibratedKappas; // paramVectorInitial;

		maturityGrid = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		LocalDate[] dates = Arrays.stream(maturityGrid).mapToObj(ttm -> {
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

		double[] forwards = Arrays.stream(maturityGrid).map(x -> equityForwardStructure.getForward(x)).toArray();

		// Initialize volatilityPoints
		strikeMatPairs = new ArrayList<>();

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

		// ACF
		acfAtLags = new ArrayList<>();
		acfAtLags.add(new Pair<>(0., 1.));
		acfAtLags.add(new Pair<>(10.,0.782428986));
		acfAtLags.add(new Pair<>(20.,0.626013798));
		acfAtLags.add(new Pair<>(30.,0.52663777));
		acfAtLags.add(new Pair<>(40.,0.474562316));
		acfAtLags.add(new Pair<>(50.,0.438915541));
		acfAtLags.add(new Pair<>(60.,0.415146521));
		acfAtLags.add(new Pair<>(70.,0.383064415));
		acfAtLags.add(new Pair<>(80.,0.33402509));
	}

	public void setDAXHestonMarchSetupSIM() {
		valuationDate = LocalDate.parse("2020-03-12");
		spot0 = 1;
		selectedParamsLNSVQD = paramVectorInitialMarchCalibrated;
		selectedParamsHeston = paramVectorHestonMarch;
		selectedParamsToCalibrate = paramVectorInitialWithPreCalibratedKappas; // paramVectorInitial;

		maturityGrid = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		LocalDate[] dates = Arrays.stream(maturityGrid).mapToObj(ttm -> {
			long days = Math.round(ttm * 365);
			return valuationDate.plusDays(days);
		}).toArray(LocalDate[]::new);

		discountDates = LNSVQDUtils.createLocalDateList(new String[]{
				"2020-03-23",
				"2020-03-30",
				"2020-04-16",
				"2020-05-18",
				"2020-06-16",
				"2020-07-16",
				"2020-08-17",
				"2020-09-16",
				"2020-10-16",
				"2020-11-16",
				"2020-12-16",
				"2021-01-18",
				"2021-02-16",
				"2021-03-16",
				"2021-09-16",
				"2022-03-16",
				"2023-03-16",
				"2024-03-18"
		});

		discountFactors = new double[]{
				1.0001509700,
				1.0002493800,
				1.0004873900,
				1.0009228300,
				1.0014010000,
				1.0018781400,
				1.0023918600,
				1.0028818000,
				1.0033950900,
				1.0039425200,
				1.0044668800,
				1.0050507800,
				1.0055805300,
				1.0060629100,
				1.0093010300,
				1.0124776200,
				1.0185304000,
				1.0238196400
		};


		disountCurve = new YieldCurve("Discount curve"
				, valuationDate
				, dayCountConvention
				, discountDates
				, discountFactors);

		forwardDates = LNSVQDUtils.createLocalDateList(new String[]{
				"2020-03-20",
				"2020-04-17",
				"2020-05-15",
				"2020-06-19",
				"2020-09-18",
				"2020-12-18",
				"2021-06-18",
				"2021-12-17",
				"2022-06-17",
				"2022-12-16",
				"2023-12-15"
		});

		forwardDiscountFactors = new double[]{
				1.000413517,
				1.002016627,
				0.999862638,
				0.999768006,
				1.004269671,
				1.007155568,
				1.006625749,
				1.008821234,
				1.008656759,
				1.009848938,
				1.008563601
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

		double[] forwards = Arrays.stream(maturityGrid).map(x -> equityForwardStructure.getForward(x)).toArray();

		// Initialize volatilityPoints
		strikeMatPairs = new ArrayList<>();

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

		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<VolatilityPoint>();

		// Create and adf volatility points
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.60, 0.723512803, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.80, 0.65872135, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.00, 0.613484389, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.20, 0.583073626, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.40, 0.564207765, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.60, 0.723512803, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.80, 0.65872135, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.00, 0.613484389, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.20, 0.583073626, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.40, 0.564207765, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.60, 0.723512803, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.80, 0.65872135, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.00, 0.613484389, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.20, 0.583073626, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.40, 0.564207765, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.60, 0.655702448, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.80, 0.597767129, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.00, 0.557298712, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.20, 0.530002894, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.40, 0.512931818, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.60, 0.593709634, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.80, 0.542123144, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.00, 0.50606706, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.20, 0.48164473, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.40, 0.466216855, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.60, 0.548502278, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.80, 0.501609555, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.00, 0.468812555, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.20, 0.446505728, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.40, 0.43227795, spot0));

		// Create volatility surface
		volatilityPointsSurface = new VolatilityPointsSurface(volatilityPoints, valuationDate, dayCountConvention);

		// ACF
		acfAtLags = new ArrayList<>();
		acfAtLags.add(new Pair<>(0., 1.));
		acfAtLags.add(new Pair<>(10.,0.782428986));
		acfAtLags.add(new Pair<>(20.,0.626013798));
		acfAtLags.add(new Pair<>(30.,0.52663777));
		acfAtLags.add(new Pair<>(40.,0.474562316));
		acfAtLags.add(new Pair<>(50.,0.438915541));
		acfAtLags.add(new Pair<>(60.,0.415146521));
		acfAtLags.add(new Pair<>(70.,0.383064415));
		acfAtLags.add(new Pair<>(80.,0.33402509));
	}

	public void setDAXHestonFebruarySetupSIM() {
		valuationDate = LocalDate.parse("2025-02-28");
		spot0 = 1;
		selectedParamsLNSVQD = paramVectorInitialFebruaryCalibrated;
		selectedParamsHeston = paramVectorHestonFebruary;
		selectedParamsToCalibrate = paramVectorInitialWithPreCalibratedKappas; // paramVectorInitial;

		maturityGrid = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		LocalDate[] dates = Arrays.stream(maturityGrid).mapToObj(ttm -> {
			long days = Math.round(ttm * 365);
			return valuationDate.plusDays(days);
		}).toArray(LocalDate[]::new);

		discountDates = LNSVQDUtils.createLocalDateList(new String[]{
				"2025-03-18",
				"2025-03-25",
				"2025-04-11",
				"2025-05-13",
				"2025-06-11",
				"2025-07-11",
				"2025-08-12",
				"2025-09-11",
				"2025-10-11",
				"2025-11-11",
				"2025-12-11",
				"2026-01-13",
				"2026-02-11",
				"2026-03-11",
				"2026-09-11",
				"2027-03-11",
				"2028-03-11",
				"2029-03-13"
		});

		double[] discountFactors = new double[]{
				0.9995301900,
				0.9990607400,
				0.9979213000,
				0.9958492400,
				0.9940334500,
				0.9922711800,
				0.9904080100,
				0.9887031900,
				0.9870389200,
				0.9853493900,
				0.9837016700,
				0.9819221300,
				0.9803808900,
				0.9789481300,
				0.9689034000,
				0.9589405200,
				0.9369840700,
				0.9136385100
		};

		disountCurve = new YieldCurve("Discount curve"
				, valuationDate
				, dayCountConvention
				, discountDates
				, discountFactors);

		forwardDates = LNSVQDUtils.createLocalDateList(new String[]{
				"2025-03-22",
				"2025-04-18",
				"2025-05-17",
				"2025-06-21",
				"2025-09-20",
				"2025-12-20",
				"2026-06-20",
				"2026-12-19",
				"2027-06-19",
				"2027-12-18",
				"2028-12-16"
		});

		double[] forwardDiscountFactors = new double[]{
				0.999056362,
				0.996571729,
				0.991297063,
				0.988033632,
				0.98145907,
				0.975199207,
				0.959150058,
				0.947124043,
				0.930399542,
				0.917067578,
				0.886655147
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

		double[] forwards = Arrays.stream(maturityGrid).map(x -> equityForwardStructure.getForward(x)).toArray();

		// Initialize volatilityPoints
		strikeMatPairs = new ArrayList<>();

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

		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<VolatilityPoint>();

		// Create and adf volatility points
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.60, 0.415786678, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.80, 0.27187155, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.00, 0.157476346, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.20, 0.150327897, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.40, 0.205672117, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.60, 0.347514351, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.80, 0.240476016, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.00, 0.155129551, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.20, 0.13254467, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.40, 0.163996755, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.60, 0.315752556, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.80, 0.226429031, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.00, 0.155941504, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.20, 0.129920126, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.40, 0.149590733, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.60, 0.295667661, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.80, 0.217692305, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.00, 0.156477694, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.20, 0.128872013, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.40, 0.140370308, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.60, 0.282470149, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.80, 0.212047644, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.00, 0.156764181, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.20, 0.128190999, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.40, 0.134168029, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.60, 0.273585698, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.80, 0.208241225, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.00, 0.157393387, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.20, 0.129904105, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.40, 0.133219151, spot0));

		// Create volatility surface
		volatilityPointsSurface = new VolatilityPointsSurface(volatilityPoints, valuationDate, dayCountConvention);

		// ACF
		acfAtLags = new ArrayList<>();
		acfAtLags.add(new Pair<>(0., 1.));
		acfAtLags.add(new Pair<>(10.,0.782428986));
		acfAtLags.add(new Pair<>(20.,0.626013798));
		acfAtLags.add(new Pair<>(30.,0.52663777));
		acfAtLags.add(new Pair<>(40.,0.474562316));
		acfAtLags.add(new Pair<>(50.,0.438915541));
		acfAtLags.add(new Pair<>(60.,0.415146521));
		acfAtLags.add(new Pair<>(70.,0.383064415));
		acfAtLags.add(new Pair<>(80.,0.33402509));
	}

	public void setBTCSetupSIM() {
		valuationDate = LocalDate.parse("2024-09-30");
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
	}

	static VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * equityForwardStructure.getForward(maturity);
		return new VolatilityPoint(maturity, strike, volatility);
	}

}
