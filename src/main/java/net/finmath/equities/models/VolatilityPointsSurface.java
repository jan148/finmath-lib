package net.finmath.equities.models;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.time.daycount.DayCountConvention;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * This class implements the volatility interfaces for a non-flat volatility surface.
 *
 * @author Andreas Grotz, Jan Berger
 */

public class VolatilityPointsSurface implements VolatilitySurface, ShiftedVolatilitySurface {
	// TODO: Check that volatility points are sorted by ExpiryDate and Strike
	private final ArrayList<VolatilityPoint> volatilityPoints;
	private final ArrayList<LocalDate> volatilityDates = new ArrayList<>();
	private final ArrayList<Double> strikes = new ArrayList<>();
	// Prices
	private final ArrayList<Double> prices = new ArrayList<>();
	private final DayCountConvention dayCountConvention;
	private final LocalDate today;
	private final double volShift;

	public VolatilityPointsSurface(ArrayList<VolatilityPoint> volatilityPoints, LocalDate today, DayCountConvention dayCountConvention, double volShift) {
		this.volatilityPoints = volatilityPoints;
		sortVolaPoints();
		volatilityPoints.forEach(volatilityPoint -> {
			volatilityDates.add(volatilityPoint.getDate());
			strikes.add(volatilityPoint.getStrike());
		});
		this.dayCountConvention = dayCountConvention;
		this.today = today;
		this.volShift = volShift;
	}

	public VolatilityPointsSurface(ArrayList<VolatilityPoint> volatilityPoints, LocalDate today, DayCountConvention dayCountConvention) {
		this(volatilityPoints, today, dayCountConvention, 0.0);
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

	public DayCountConvention getDayCountConvention() {
		return this.dayCountConvention;
	}


	public LocalDate getToday() {
		return this.today;
	}

	public double getShift() {
		return volShift;
	}

	@Override
	public ShiftedVolatilitySurface getShiftedSurface(double shift) {
		assert volShift == 0.0 : "Surface is already shifted";
		return new VolatilityPointsSurface(volatilityPoints, today, dayCountConvention, shift);
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

	public void printVolSurfaceForOutput() {
		System.out.println("Date \t Strike \t Volatility");
		for(VolatilityPoint volatilityPoint : volatilityPoints) {
			System.out.println(volatilityPoint.getDate() + "\t" + volatilityPoint.getStrike() + "\t" + volatilityPoint.getVolatility());
		}
	}

	private void sortVolaPoints() {
		Comparator<VolatilityPoint> comparator = new Comparator<VolatilityPoint>() {
			@Override
			public int compare(VolatilityPoint o1, VolatilityPoint o2) {
				if(o1.getDate().isBefore(o2.getDate())) {
					return -1;
				} else if(o1.getDate().isEqual(o2.getDate()) && o1.getStrike() < o2.getStrike()) {
					return -1;
				} else {
					return 1;
				}
			}
		};
		Collections.sort(volatilityPoints, comparator);
	}

}
