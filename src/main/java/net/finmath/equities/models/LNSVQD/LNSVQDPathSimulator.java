package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.time.LocalDate;

public abstract class LNSVQDPathSimulator extends PathSimulator{
	LNSVQDModel lnsvqdModel;
	boolean isBackwardEuler = false;
	// y = e^x
	UnivariateFunction zeta = y -> 1 / y * lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() - y * lnsvqdModel.getKappa2()
			- lnsvqdModel.getKappa1() + lnsvqdModel.getKappa2() * lnsvqdModel.getTheta() - 0.5 * lnsvqdModel.getTotalInstVar();

	public LNSVQDPathSimulator(LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure, int numberOfPaths, double[] timeGrid, double[] maturities, LNSVQDModel lnsvqdModel, Boolean isBackwardEuler) {
		this.spotDate = spotDate;
		this.discountCurve = discountCurve;
		this.equityForwardStructure = equityForwardStructure;
		this.numberOfPaths = numberOfPaths;
		this.timeGrid = timeGrid;
		this.maturities = maturities;
		this.lnsvqdModel = lnsvqdModel;
		this.isBackwardEuler = isBackwardEuler;
	}
}
