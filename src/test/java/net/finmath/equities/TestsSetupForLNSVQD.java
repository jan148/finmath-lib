package net.finmath.equities;

import net.finmath.equities.marketdata.AffineDividend;
import net.finmath.equities.marketdata.AffineDividendStream;
import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.BuehlerDividendForwardStructure;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQD.LNSVQDCallPriceSimulator;
import net.finmath.equities.models.LNSVQD.LNSVQDModelAnalyticalPricer;
import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.equities.models.VolatilityPointsSurface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.time.daycount.DayCountConvention;
import net.finmath.time.daycount.DayCountConvention_ACT_365;
import org.junit.jupiter.api.BeforeEach;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Math.E;

public abstract class TestsSetupForLNSVQD {
	/**
	 * Time params
	 */
	LocalDate valuationDate = LocalDate.parse("2024-09-30");
	double spot0 = 1; // 19324.93;

	/**
	 * Create analytical model
	 */
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

	double[] selectedParams = paramVectorCalibrated;

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

	LNSVQDModelAnalyticalPricer lnsvqdModelAnalyticalPricer =
			new LNSVQDModelAnalyticalPricer(
					spot0
					, selectedParams[0]
					, selectedParams[1]
					, selectedParams[2]
					, selectedParams[3]
					, selectedParams[4]
					, selectedParams[5]
					, 0, valuationDate, equityForwardStructure);

	/**
	 * Create simulation model (not finmath)
	 */
	int numberOfPaths = 100000;
	// TODO: Discounts dates should be decoupled from maturities
	/*double[] maturityGrid = Arrays.stream(discountDates)
			.mapToDouble(date -> dayCountConvention.getDaycountFraction(valuationDate, date))
			.toArray();*/
	double[] maturityGrid = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
	int numberPointsToInsert = (int) (maturityGrid[maturityGrid.length - 1] * 365 * 2 - maturityGrid.length);
	List<Double> timeGridForSimulationList;

	{
		try {
			timeGridForSimulationList = LNSVQDUtils.addTimePointsToArray
					(maturityGrid, numberPointsToInsert, 0, maturityGrid[maturityGrid.length - 1] - 1E-18, true);
		} catch(Exception e) {
			throw new RuntimeException(e);
		}
	}

	double[] timeGridForSimulation = timeGridForSimulationList.stream()
			.mapToDouble(Double::doubleValue)
			.toArray();
	LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator =
			new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGridForSimulation, false);

	/**
	 * Create simulation model (finmath)
	 */
	TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(timeGridForSimulation);
	RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

	/**
	 * Declare volatility surface; will be insatntiated later
	 */
	VolatilityPointsSurface volatilityPointsSurface;

	//-----------------------------------------------------------------------------------------------------------------

	/**
	 * ***************************************************+
	 * Init method
	 * ***************************************************+
	 */
	@BeforeEach
	public void init() {
		setTargetSurfaceHeston();
	}

	//-----------------------------------------------------------------------------------------------------------------

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


	public void setTargetSurfaceHeston() {
		double[] ttms = new double[]{0.25, 0.5, 0.75, 1, 1.25, 1.5};
		LocalDate[] dates = Arrays.stream(ttms).mapToObj(ttm -> {
			long days = Math.round(ttm * 365);
			return valuationDate.plusDays(days);
		}).toArray(LocalDate[]::new);

		// Initialize volatilityPoints
		ArrayList<VolatilityPoint> volatilityPoints = new ArrayList<>();

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

		// Redefine the maturityGrid and the simulation grid
		maturityGrid = Arrays.stream(dates).mapToDouble(date -> dayCountConvention.getDaycountFraction(valuationDate, date)).toArray();
		int numberPointsToInsert = (int) (maturityGrid[maturityGrid.length - 1] * 365 * 2 - maturityGrid.length);
		List<Double> timeGridForSimulationList;
		try {
			timeGridForSimulationList = LNSVQDUtils.addTimePointsToArray
					(maturityGrid, numberPointsToInsert, 0, maturityGrid[maturityGrid.length - 1] - 1E-18, true);
		} catch(Exception e) {
			throw new RuntimeException(e);
		}
		timeGridForSimulation = timeGridForSimulationList.stream()
				.mapToDouble(Double::doubleValue)
				.toArray();
	}

	private VolatilityPoint makeVolatilityPoint(String date, double percentage, double volatility, double spot) {
		LocalDate maturity = LocalDate.parse(date);
		double strike = percentage * spot * equityForwardStructure.getForward(maturity);
		return new VolatilityPoint(maturity, strike, volatility);
	}

}
