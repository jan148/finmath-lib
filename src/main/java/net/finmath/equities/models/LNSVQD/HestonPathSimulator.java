package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.time.LocalDate;

public abstract class HestonPathSimulator extends PathSimulator{
	/**
	 * Model parameters under the EMM
	 */
	double sigma0;
	double kappa;
	double theta;
	double epsilon;
	double rho;

	/**
	 * Market observables
	 */

	public HestonPathSimulator(LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure, int numberOfPaths, double[] timeGrid, double[] maturities, double sigma0, double kappa, double theta, double epsilon, double rho) {
		this.spotDate = spotDate;
		this.discountCurve = discountCurve;
		this.equityForwardStructure = equityForwardStructure;
		this.numberOfPaths = numberOfPaths;
		this.timeGrid = timeGrid;
		this.maturities = maturities;
		this.sigma0 = sigma0;
		this.kappa = kappa;
		this.theta = theta;
		this.epsilon = epsilon;
		this.rho = rho;
	}

}
