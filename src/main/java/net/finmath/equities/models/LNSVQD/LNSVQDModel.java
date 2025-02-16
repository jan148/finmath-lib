package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.AffineDividendStream;
import net.finmath.equities.marketdata.FlatYieldCurve;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.model.AbstractProcessModel;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.daycount.DayCountConvention;
import net.finmath.time.daycount.DayCountConvention_ACT_365;

import java.time.LocalDate;
import java.util.Map;
import java.util.function.Function;

/**
 * TODO:
 * 1. Delete the quantitites that are derived from the parameter values create methods for them instead (less error-prone)
 */

public class LNSVQDModel extends AbstractProcessModel {
	/**
	 * Model parameters under the EMM
	 */
	protected final double spot0;
	double sigma0;
	double kappa1;
	double kappa2;
	double theta;
	double beta;
	double epsilon;
	double totalInstVar;

	/**
	 * Market observables
	 */
	LocalDate spotDate;
	EquityForwardStructure equityForwardStructure;

	/**
	 * Transformed inital values
	 */
	protected double X0, Y0, I0;

	/**
	 * Random variable factory
	 */
	protected final RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();
	protected static final RandomVariable ZERO = new Scalar(0.0);

	/**
	 * Functions used for the discretization scheme of the LNSVQD model
	 */
	public Function<RandomVariable, RandomVariable> zeta = new Function<RandomVariable, RandomVariable>() {
		@Override
		public RandomVariable apply(RandomVariable randomVariable) {
			return randomVariable.mult(-1).exp().mult(getKappa1() * getTheta()).sub(randomVariable.exp().mult(getKappa2()))
					.add(-getKappa1() + getKappa2() * getTheta() - 0.5 * getTotalInstVar());
		}

	};

	public Function<RandomVariable, RandomVariable> zeta1stDerivative = new Function<RandomVariable, RandomVariable>() {
		@Override
		public RandomVariable apply(RandomVariable randomVariable) {
			return randomVariable.mult(-1).exp().mult(getKappa1() * getTheta()).mult(-1)
					.sub(randomVariable.exp().mult(getKappa2()));
		}
	};

	public LNSVQDModel(double spot0, double sigma0, double kappa1, double kappa2, double theta, double beta, double epsilon, double I0, LocalDate spotDate, EquityForwardStructure equityForwardStructure) {
		super();

		// Perform necessary checks
		checkMartingalityOfDiscountedAssetProcess(kappa2, beta);

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

		// EFS
		this.spotDate = spotDate;
		this.equityForwardStructure = equityForwardStructure;
	}

	public LNSVQDModel(double spot0, double sigma0, double kappa1, double kappa2, double theta, double beta, double epsilon, double I0, LocalDate spotDate) {
		super();

		// Perform necessary checks
		checkMartingalityOfDiscountedAssetProcess(kappa2, beta);

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

		// EFS; create a EFS with a flat yield curve and 365 dayCountConvention
		this.spotDate = spotDate;
		this.equityForwardStructure = new EquityForwardStructure() {
			@Override
			public DividendModelType getDividendModel() {
				return null;
			}

			@Override
			public LocalDate getValuationDate() {
				return spotDate;
			}

			@Override
			public double getSpot() {
				return 0;
			}

			@Override
			public YieldCurve getRepoCurve() {
				return new FlatYieldCurve(getValuationDate(), 0.00, new DayCountConvention_ACT_365());
			}

			@Override
			public AffineDividendStream getDividendStream() {
				return null;
			}

			@Override
			public EquityForwardStructure cloneWithNewSpot(double newSpot) {
				return null;
			}

			@Override
			public EquityForwardStructure cloneWithNewDate(LocalDate newDate) {
				return null;
			}

			@Override
			public double getGrowthDiscountFactor(double startTime, double endTime) {
				return 0;
			}

			@Override
			public double getGrowthDiscountFactor(LocalDate startDate, LocalDate endDate) {
				return 0;
			}

			@Override
			public double getFutureDividendFactor(double valTime) {
				return 0;
			}

			@Override
			public double getFutureDividendFactor(LocalDate valDate) {
				return 0;
			}

			@Override
			public double getForward(double expiryTime) {
				return 0;
			}

			@Override
			public double getForward(LocalDate expiryDate) {
				return 0;
			}

			@Override
			public double getDividendAdjustedStrike(double strike, double expiryTime) {
				return 0;
			}

			@Override
			public double getDividendAdjustedStrike(double strike, LocalDate expiryDate) {
				return 0;
			}

			@Override
			public double getLogMoneyness(double strike, double expiryTime) {
				return 0;
			}

			@Override
			public double getLogMoneyness(double strike, LocalDate expiryDate) {
				return 0;
			}
		};
	}

	/**
	 * ***************************************************+
	 * Getters
	 * ***************************************************+
	 */
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

	public double getX0() {
		return X0;
	}

	public double getY0() {
		return Y0;
	}

	public double getI0() {
		return I0;
	}

	public double getTotalInstVar() {
		return totalInstVar;
	}

	public double getRiskFreeRate(double time) {
		return equityForwardStructure.getRepoCurve().getRate(time);
	}

	/**
	 * ***************************************************+
	 * Setters
	 * ***************************************************+
	 */
	public void setVolatilityParameters(double[] parameterVector) {
		// Perform necessary checks
		// checkMartingalityOfDiscountedAssetProcess(kappa2, beta);

		this.sigma0 = parameterVector[0];
		this.kappa1 = parameterVector[1];
		this.kappa2 = parameterVector[2];
		this.theta = parameterVector[3];
		this.beta = parameterVector[4];
		this.epsilon = parameterVector[5];
		this.totalInstVar = beta * beta + epsilon * epsilon;

		this.X0 = Math.log(this.spot0);
		this.Y0 = sigma0 - theta;
		this.I0 = I0;
	}

	/**
	 * ***************************************************+
	 * SECTION 1: Check for parameter conditions
	 * ***************************************************+
	 */

	/**
	 * Assumption 3.1: κ1 ≥ 0, κ2 ≥ 0, θ > 0, ϑ ≥ 0.
	 */
	private void volatilityProcessParametersCondition() {
		if(!(kappa1 >= 0 && kappa2 >= 0 && theta > 0 && Math.sqrt(totalInstVar) >= 0)) {
			throw new IllegalStateException();
		}
	}

	/**
	 * Condition 1 in Theorem 3.7: κ2 ≥ beta.
 	 */
	private void checkMartingalityOfDiscountedAssetProcess(double kappa2, double beta) {
		if(kappa2 < beta) {
			// throw new IllegalStateException("κ2 < beta. Martingale condition violated!");
		}
	}

	/**
	 * ***************************************************+
	 * SECTION 2: Simulation
	 * ***************************************************+
	 */

	@Override
	public int getNumberOfComponents() {
		return 2;
	}

	/**
	 * Map from (S, sigma) to (ln(S / M), ln sigma)
	 */
	@Override
	public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		double time = process.getTime(timeIndex);
		if(componentIndex == 0) {
			return randomVariable.div(getNumeraire(null, time)).log();
		}
		else if(componentIndex == 1) {
			return randomVariable.log();
		}
		else {
			throw new UnsupportedOperationException("Component " + componentIndex + " does not exist.");
		}
	}

	/**
	 * Map from (ln(S / M), ln sigma) to (S, sigma)
	 */
	@Override
	public RandomVariable applyStateSpaceTransformInverse(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		double time = process.getTime(timeIndex);
		if(componentIndex == 0) {
			return randomVariable.exp().mult(getNumeraire(null, time));
		}
		else if(componentIndex == 1) {
			return randomVariable.exp();
		}
		else {
			throw new UnsupportedOperationException("Component " + componentIndex + " does not exist.");
		}
	}

	@Override
	public RandomVariable[] getInitialState(MonteCarloProcess process) {
		RandomVariable[] initialValueVector = new RandomVariable[2];
		initialValueVector[0] = randomVariableFactory.createRandomVariable(spot0);
		initialValueVector[1] = randomVariableFactory.createRandomVariable(sigma0);
		return initialValueVector;
	}

	@Override
	public RandomVariable getNumeraire(MonteCarloProcess process, double time) {
		return getRandomVariableForConstant(1. / equityForwardStructure.getRepoCurve().getDiscountFactor(time));
	}

	/**
	 * @param realizationAtTimeIndex: Realizations of (S, \sigma), i.e. of the original process
	 * @return The drift of the transformed process, i.e. the one used for the MC-discretization scheme;
	 *
	 * NOTE: We only define the scheme for the asset process; The implicit scheme for the volatility process is defined
	 * in the class repsonsible for the discretization scheme
	 */
	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		RandomVariable stochasticVolatility = realizationAtTimeIndex[1];
		// RandomVariable lnStochasticVolatility = applyStateSpaceTransform(process, timeIndex, 1, realizationAtTimeIndex[1]);
		// TODO: Check the following formulas
		RandomVariable driftAsset = stochasticVolatility.pow(2).mult(-0.5);
		// RandomVariable driftVolatility = zeta.apply(lnStochasticVolatility);

		return new RandomVariable[]{driftAsset, ZERO};
	}

	@Override
	public int getNumberOfFactors() {
		return 2;
	}

	/**
	 * @return The factor loadings of the transformed process, i.e. the one used for the MC-discretization scheme
	 */
	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		RandomVariable stochasticVolatility = realizationAtTimeIndex[1];
		final RandomVariable[] factorLoadings = new RandomVariable[2];
		if(componentIndex == 0) {
			factorLoadings[0] = stochasticVolatility;
			factorLoadings[1] = ZERO;
		}
		else if(componentIndex == 1) {
			factorLoadings[0] = getRandomVariableForConstant(getBeta());
			factorLoadings[1] = getRandomVariableForConstant(getEpsilon());
		}
		else {
			throw new UnsupportedOperationException("Component " + componentIndex + " does not exist.");
		}
		// Return factor loadings
		return factorLoadings;
	}

	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return randomVariableFactory.createRandomVariable(value);
	}

	public RandomVariable getRandomVariableForArray(double[] values) {
		return randomVariableFactory.createRandomVariable(-1, values);
	}

	// TODO
	@Override
	public ProcessModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		return null;
	}

}
