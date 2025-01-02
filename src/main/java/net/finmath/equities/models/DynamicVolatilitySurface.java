package net.finmath.equities.models;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.functions.AnalyticFormulas;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;

/**
 * This class implements the volatility interfaces for a non-flat volatility surface.
 *
 * @author Andreas Grotz, Jan Berger
 */

public class DynamicVolatilitySurface implements VolatilitySurface, ShiftedVolatilitySurface {
	// TODO: Check that volatility points are sorted by ExpiryDate and Strike
	private final ArrayList<VolatilityPoint> volatilityPoints;
	private final ArrayList<LocalDate> volatilityDates = new ArrayList<>();
	private final ArrayList<Double> strikes = new ArrayList<>();
	// Prices
	private final ArrayList<Double> prices = new ArrayList<>();
	private final LocalDate today;
	private final double volShift;

	public DynamicVolatilitySurface(ArrayList<VolatilityPoint> volatilityPoints) {
		this(volatilityPoints, 0.0);
	}

	public DynamicVolatilitySurface(ArrayList<VolatilityPoint> volatilityPoints, double volShift) {
		this.volatilityPoints = volatilityPoints;
		volatilityPoints.forEach(volatilityPoint -> {
			volatilityDates.add(volatilityPoint.getDate());
			strikes.add(volatilityPoint.getVolatility());
		});
		this.today = volatilityDates.get(0);
		this.volShift = volShift;
	}

	public int getNumberOfVolatilityPoints() {
		return volatilityPoints.size();
	}
	public ArrayList<VolatilityPoint> getVolatilityPoints() {
		return volatilityPoints;
	}

	public ArrayList<LocalDate> getVolatilityDates() {
		return volatilityDates;
	}

	public ArrayList<Double> getStrikes() {
		return strikes;
	}

	public LocalDate getToday() {
		return this.today;
	}

	public double getShift() {
		return volShift;
	}

	public ArrayList<Double> getPrices(double initalStockValue, double riskFreeRate) {
		ArrayList<Double> prices = new ArrayList<>();

		for(VolatilityPoint volatilityPoint : volatilityPoints) {
			LocalDate date = volatilityPoint.getDate();
			long period = ChronoUnit.DAYS.between(today, date);
			double ttm = period / 365.; // TODO: Get the correct denominator

			// Fetch point from the surface
			double strike = volatilityPoint.getStrike();
			double empiricalVolatility = volatilityPoint.getVolatility();
			double price = AnalyticFormulas.blackScholesOptionValue(initalStockValue, riskFreeRate, empiricalVolatility, ttm, strike);
			prices.add(price);
		}

		return prices;
	}

	@Override
	public ShiftedVolatilitySurface getShiftedSurface(double shift) {
		assert volShift == 0.0 : "Surface is already shifted";
		return new DynamicVolatilitySurface(volatilityPoints, shift);
	}

	@Override
	public double getVolatility(
			double strike,
			LocalDate expiryDate,
			EquityForwardStructure currentForwardStructure) {
		int index = volatilityDates.indexOf(expiryDate);
		ArrayList<Double> subList = new ArrayList<>(strikes.subList(index, strikes.size()));
		index += subList.indexOf(strike);
		return volatilityPoints.get(index).getVolatility() + volShift;
	}

	// TODO
	@Override
	public double getVolatility(double strike, double timeToMaturity, EquityForwardStructure currentForwardStructure) {
		return 0;
	}

	// TODO
	@Override
	public double getLocalVolatility(
			double strike,
			LocalDate expiryDate,
			EquityForwardStructure currentForwardStructure,
			double strikeShift,
			double timeShift) {
		return 0;
	}

	// TODO
	@Override
	public double getLocalVolatility(
			double logStrike,
			double timeToMaturity,
			EquityForwardStructure currentForwardStructure,
			double strikeShift,
			double timeShift) {
		return 0;
	}

	// TODO
	@Override
	public void calibrate(EquityForwardStructure forwardStructure, ArrayList<VolatilityPoint> volaPoints) {}

}
