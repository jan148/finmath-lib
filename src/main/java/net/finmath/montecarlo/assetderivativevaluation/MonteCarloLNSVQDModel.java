package net.finmath.montecarlo.assetderivativevaluation;

import net.finmath.equities.models.LNSVQDModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.assetderivativevaluation.models.MertonModel;
import net.finmath.montecarlo.process.LNSVQDDiscretizationScheme;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.Map;

/**
 * This class glues together a <code>MertonModel</code> and a Monte-Carlo implementation of a <code>MonteCarloProcessFromProcessModel</code>, namely <code>EulerSchemeFromProcessModel</code>,
 * and forms a Monte-Carlo implementation of the Merton model by implementing <code>AssetModelMonteCarloSimulationModel</code>.
 *
 * The model is
 * \[
 * 	dS = \mu S dt + \sigma S dW + S dJ, \quad S(0) = S_{0},
 * \]
 * \[
 * 	dN = r N dt, \quad N(0) = N_{0},
 * \]
 * where \( W \) is Brownian motion and \( J \)  is a jump process (compound Poisson process).
 *
 * The process \( J \) is given by \( J(t) = \sum_{i=1}^{N(t)} (Y_{i}-1) \), where
 * \( \log(Y_{i}) \) are i.i.d. normals with mean \( a - \frac{1}{2} b^{2} \) and standard deviation \( b \).
 * Here \( a \) is the jump size mean and \( b \) is the jump size std. dev.
 *
 * For details on the construction of the model see {@link MertonModel}.
 *
 * @author Christian Fries, Jan Berger
 * @see MertonModel
 * @see MonteCarloProcess The interface for numerical schemes.
 * @see net.finmath.montecarlo.model.ProcessModel The interface for models provinding parameters to numerical schemes.
 * @version 1.0
 */
public class MonteCarloLNSVQDModel implements AssetModelMonteCarloSimulationModel {

	private final LNSVQDModel model;
	private final MonteCarloProcess process;

	private final double initialValue;
	private final int seed;

	/**
	 * Create a Monte-Carlo simulation using given time discretization and given parameters.
	 *
	 * @param seed The seed used for the random number generator.
	 */
	public MonteCarloLNSVQDModel(
			final LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme,
			final int seed
			) {
		super();

		this.initialValue = lnsvqdDiscretizationScheme.getModel().getSpot0();
		this.seed = seed;

		// Create the model
		model = lnsvqdDiscretizationScheme.getModel();

		// Create a corresponding MC process
		process = lnsvqdDiscretizationScheme;
	}

	@Override
	public LocalDateTime getReferenceDate() {
		throw new UnsupportedOperationException("This model does not provide a reference date. Reference dates will be mandatory in a future version.");
	}

	@Override
	public RandomVariable getAssetValue(final double time, final int assetIndex) throws CalculationException {
		final int timeIndex = getTimeIndex(time);
		if(timeIndex < 0) {
			throw new IllegalArgumentException("The model does not provide an interpolation of simulation time (time given was " + time + ").");
		}

		return getAssetValue(timeIndex, assetIndex);
	}

	@Override
	public RandomVariable getAssetValue(final int timeIndex, final int assetIndex) throws CalculationException {
		return process.getProcessValue(timeIndex, assetIndex);
	}

	@Override
	public RandomVariable getNumeraire(final int timeIndex) throws CalculationException {
		final double time = getTime(timeIndex);

		return model.getNumeraire(process, time);
	}

	@Override
	public RandomVariable getNumeraire(final double time) throws CalculationException {
		return model.getNumeraire(process, time);
	}

	@Override
	public RandomVariable getMonteCarloWeights(final double time) throws CalculationException {
		return getMonteCarloWeights(getTimeIndex(time));
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationModel#getNumberOfAssets()
	 */
	@Override
	public int getNumberOfAssets() {
		return 1;
	}

	// TODO: Finish the method
	@Override
	public AssetModelMonteCarloSimulationModel getCloneWithModifiedData(final Map<String, Object> dataModified) {
		/*
		 * Determine the new model parameters from the provided parameter map.
		 */
		/*final double	newInitialTime	= dataModified.get("initialTime") != null	? ((Number)dataModified.get("initialTime")).doubleValue() : getTime(0);
		final double	newInitialValue	= dataModified.get("initialValue") != null	? ((Number)dataModified.get("initialValue")).doubleValue() : initialValue;
		final double	newRiskFreeRate	= dataModified.get("riskFreeRate") != null	? ((Number)dataModified.get("riskFreeRate")).doubleValue() : model.getRiskFreeRate().doubleValue();
		final double	newVolatility	= dataModified.get("volatility") != null	? ((Number)dataModified.get("volatility")).doubleValue()	: model.getVolatility().doubleValue();
		final double	newJumpIntensity	= dataModified.get("jumpIntensity") != null	? ((Number)dataModified.get("jumpIntensity")).doubleValue()	: model.getJumpIntensity().doubleValue();
		final double	newJumpSizeMean		= dataModified.get("jumpSizeMean") != null	? ((Number)dataModified.get("jumpSizeMean")).doubleValue()	: model.getVolatility().doubleValue();
		final double	newJumpSizeStdDev	= dataModified.get("jumpSizeStdDev") != null	? ((Number)dataModified.get("jumpSizeStdDev")).doubleValue()	: model.getVolatility().doubleValue();
		final int		newSeed				= dataModified.get("seed") != null			? ((Number)dataModified.get("seed")).intValue()				: seed;

		return new MonteCarloLNSVQDModel(process.getTimeDiscretization().getTimeShiftedTimeDiscretization(newInitialTime-getTime(0)), process.getNumberOfPaths(), newSeed, newInitialValue, newRiskFreeRate, newVolatility, newJumpIntensity, newJumpSizeMean, newJumpSizeStdDev);*/
		return null;

	}

	@Override
	public AssetModelMonteCarloSimulationModel getCloneWithModifiedSeed(final int seed) {
		final Map<String, Object> dataModified = new HashMap<>();
		dataModified.put("seed", new Integer(seed));
		return getCloneWithModifiedData(dataModified);
	}

	@Override
	public int getNumberOfPaths() {
		return process.getNumberOfPaths();
	}

	@Override
	public TimeDiscretization getTimeDiscretization() {
		return process.getTimeDiscretization();
	}

	@Override
	public double getTime(final int timeIndex) {
		return process.getTime(timeIndex);
	}

	@Override
	public int getTimeIndex(final double time) {
		return process.getTimeIndex(time);
	}

	@Override
	public RandomVariable getRandomVariableForConstant(final double value) {
		return model.getRandomVariableForConstant(value);
	}

	@Override
	public RandomVariable getMonteCarloWeights(final int timeIndex) throws CalculationException {
		return process.getMonteCarloWeights(timeIndex);
	}
}
