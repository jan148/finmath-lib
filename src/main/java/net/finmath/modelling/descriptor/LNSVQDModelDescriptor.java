/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 09.02.2018
 */

package net.finmath.modelling.descriptor;

import net.finmath.equities.marketdata.AffineDividendStream;
import net.finmath.equities.marketdata.FlatYieldCurve;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.Black76Model;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.daycount.DayCountConvention_ACT_365;

import java.time.LocalDate;
import java.util.Map;
import java.util.function.Function;

/**
 * Model descriptor for the log-normal stcohastic volatility model with quadratic drift,
 * see https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2522425.
 *
 * @author Jan Berger
 *
 * @version 1.0
 */
public class LNSVQDModelDescriptor implements AssetModelDescriptor {
	/**
	 * Model parameters under the EMM
	 */
	protected final double spot0;
	protected final double sigma0;
	protected final double kappa1;
	protected final double kappa2;
	protected final double theta;
	protected final double beta;
	protected final double epsilon;
	protected final double totalInstVar;

	/**
	 * Transformed initial values
	 */
	protected final double X0, Y0, I0;

	/**
	 * Market observables
	 */
	protected final LocalDate spotDate;
	protected final YieldCurve discountCurve;
	protected final EquityForwardStructure equityForwardStructure;


	public LNSVQDModelDescriptor(double spot0, double sigma0, double kappa1, double kappa2, double theta, double beta
			, double epsilon, double I0, LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure) {
		super();

		this.spot0 = spot0;
		this.sigma0 = sigma0;
		this.kappa1 = kappa1;
		this.kappa2 = kappa2;
		this.theta = theta;
		this.beta = beta;
		this.epsilon = epsilon;
		this.totalInstVar = beta * beta + epsilon * epsilon;

		this.X0 = Math.log(this.spot0);
		this.Y0 = sigma0 - theta;
		this.I0 = I0;

		this.spotDate = spotDate;
		this.discountCurve=discountCurve;
		this.equityForwardStructure=equityForwardStructure;
	}

	@Override
	public Integer version() {
		return 1;
	}

	@Override
	public String name() {
		return "LNSVQD model";
	}

	public double getSpot0() {
		return spot0;
	}

	public double getSigma0() {
		return sigma0;
	}

	public double getKappa1() {
		return kappa1;
	}

	public double getKappa2() {
		return kappa2;
	}

	public double getTheta() {
		return theta;
	}

	public double getBeta() {
		return beta;
	}

	public double getEpsilon() {
		return epsilon;
	}

	public double getTotalInstVar() {
		return totalInstVar;
	}

	public double getX0() {
		return X0;
	}

	public double getY0() {
		return Y0;
	}

	public double getI0() {
		return I0;
	}

	public LocalDate getSpotDate() { return spotDate; }

	public YieldCurve getDiscountCurve() { return discountCurve;}

	public EquityForwardStructure getEquityForwardStructure() {
		return equityForwardStructure;
	}
}
