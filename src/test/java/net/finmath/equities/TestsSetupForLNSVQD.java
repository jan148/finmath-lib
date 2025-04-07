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
	static double[] selectedParamsLNSVQD;
	static double[] selectedParamsHeston;
	static double[] selectedParamsToCalibrate;

	/**
	 * Other
	 */
	static DayCountConvention dayCountConvention = new DayCountConvention_ACT_365();
	static List<Pair<Double, Double>> acfAtLags;
	static List<Pair<Double, Double>> lnSigmaSteadyStateDist;

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
	static int numberOfPaths = 1000;

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
		selectedParamsLNSVQD = ParametersLNSVQD.paramInitDAXSepCalib;
		selectedParamsHeston = ParametersLNSVQD.paramVectorHestonDAXSep;
		selectedParamsToCalibrate = ParametersLNSVQD.paramInitDAXSep;

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

		// Steady state distribution
		lnSigmaSteadyStateDist = new ArrayList<>();
		lnSigmaSteadyStateDist.add(new Pair<>(-2.2, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-2.1, 0.017767478));
		lnSigmaSteadyStateDist.add(new Pair<>(-2., 0.076477404));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.9, 0.129393588));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.8, 0.1201236));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.7, 0.129007339));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.6, 0.123599846));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.5, 0.09578988));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.4, 0.100038625));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.3, 0.065276168));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.2, 0.059868675));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.1, 0.0397837));
		lnSigmaSteadyStateDist.add(new Pair<>(-1., 0.01622248));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.9, 0.010428737));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.8, 0.005021244));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.7, 0.002703747));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.6, 0.002703747));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.5, 0.001544998));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.4, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.3, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.2, 0.001158749));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.1, 0.00077249 ));
	}

	public void setDAXHestonMarchSetupSIM() {
		valuationDate = LocalDate.parse("2020-03-12");
		spot0 = 9254.07;
		selectedParamsLNSVQD = ParametersLNSVQD.paramInitDAXMarchCalib;
		selectedParamsHeston = ParametersLNSVQD.paramVectorHestonDAXMarch;
		selectedParamsToCalibrate = ParametersLNSVQD.paramInitDAXMarch;

		maturityGrid = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		LocalDate[] dates = Arrays.stream(maturityGrid).mapToObj(ttm -> {
			long days = Math.round(ttm * 365);
			return valuationDate.plusDays(days);
		}).toArray(LocalDate[]::new);

		discountDates = LNSVQDUtils.createLocalDateList(new String[]{
				"2020-03-19",
				"2020-03-26",
				"2020-04-14",
				"2020-05-12",
				"2020-06-12",
				"2020-07-13",
				"2020-08-12",
				"2020-09-14",
				"2020-10-12",
				"2020-11-12",
				"2020-12-14",
				"2021-01-12",
				"2021-02-12",
				"2021-03-12",
				"2021-09-13",
				"2022-03-14",
				"2023-03-13",
				"2024-03-12"
		});

		discountFactors = new double[]{
				1.0001114324,
				1.0002151063,
				1.0005360772,
				1.0010092075,
				1.0015382526,
				1.0021125334,
				1.0026675971,
				1.0032760273,
				1.0038369459,
				1.0044707285,
				1.0051353471,
				1.0057276188,
				1.0063636702,
				1.0069114705,
				1.0107826153,
				1.0145816867,
				1.0213373294,
				1.0274204190
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
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.60, 0.785448623, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 0.80, 0.64250993, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.00, 0.522559669, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.20, 0.4330327, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[0].toString(), 1.40, 0.389166585, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.60, 0.63284719, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 0.80, 0.525241632, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.00, 0.435594447, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.20, 0.365653625, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[1].toString(), 1.40, 0.322900301, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.60, 0.555076443, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 0.80, 0.464934416, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.00, 0.390344892, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.20, 0.33129643, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[2].toString(), 1.40, 0.292097563, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.60, 0.498391651, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 0.80, 0.4213057, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.00, 0.358264968, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.20, 0.308548525, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[3].toString(), 1.40, 0.274590694, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.60, 0.460228523, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 0.80, 0.392186506, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.00, 0.337122623, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.20, 0.293820364, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[4].toString(), 1.40, 0.263493011, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.60, 0.429366647, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 0.80, 0.369927047, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.00, 0.32209488, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.20, 0.284219166, spot0));
		volatilityPoints.add(makeVolatilityPoint(dates[5].toString(), 1.40, 0.25657709, spot0));

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

		// Steady state distribution
		lnSigmaSteadyStateDist = new ArrayList<>();
		lnSigmaSteadyStateDist.add(new Pair<>(-2.2, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-2.1, 0.017767478));
		lnSigmaSteadyStateDist.add(new Pair<>(-2., 0.076477404));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.9, 0.129393588));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.8, 0.1201236));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.7, 0.129007339));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.6, 0.123599846));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.5, 0.09578988));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.4, 0.100038625));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.3, 0.065276168));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.2, 0.059868675));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.1, 0.0397837));
		lnSigmaSteadyStateDist.add(new Pair<>(-1., 0.01622248));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.9, 0.010428737));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.8, 0.005021244));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.7, 0.002703747));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.6, 0.002703747));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.5, 0.001544998));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.4, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.3, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.2, 0.001158749));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.1, 0.00077249 ));
	}

	public void setDAXHestonFebruarySetupSIM() {
		valuationDate = LocalDate.parse("2025-02-28");
		spot0 = 1;
		selectedParamsLNSVQD = ParametersLNSVQD.paramInitDAXFebCalib;
		selectedParamsHeston = ParametersLNSVQD.paramVectorHestonDAXFeb;
		selectedParamsToCalibrate = ParametersLNSVQD.paramInitDAXFeb;

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

		// Steady state distribution
		lnSigmaSteadyStateDist = new ArrayList<>();
		lnSigmaSteadyStateDist.add(new Pair<>(-2.2, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-2.1, 0.017767478));
		lnSigmaSteadyStateDist.add(new Pair<>(-2., 0.076477404));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.9, 0.129393588));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.8, 0.1201236));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.7, 0.129007339));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.6, 0.123599846));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.5, 0.09578988));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.4, 0.100038625));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.3, 0.065276168));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.2, 0.059868675));
		lnSigmaSteadyStateDist.add(new Pair<>(-1.1, 0.0397837));
		lnSigmaSteadyStateDist.add(new Pair<>(-1., 0.01622248));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.9, 0.010428737));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.8, 0.005021244));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.7, 0.002703747));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.6, 0.002703747));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.5, 0.001544998));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.4, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.3, 0.000772499));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.2, 0.001158749));
		lnSigmaSteadyStateDist.add(new Pair<>(-0.1, 0.00077249 ));
	}

	public void setBTCSetupSIM() {
		valuationDate = LocalDate.parse("2024-09-30");
		spot0 = 67843.219;
		selectedParamsLNSVQD = ParametersLNSVQD.paramVectorBitcoin;
		selectedParamsHeston = ParametersLNSVQD.paramVectorHestonDAXSep;

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
