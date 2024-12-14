package net.finmath.montecarlo.process;

import net.finmath.concurrency.FutureWrapper;
import net.finmath.equities.models.LNSVQDModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.montecarlo.process.MonteCarloProcessFromProcessModel;
import net.finmath.rootfinder.NewtonsMethod;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.*;
import java.util.function.BiFunction;
import java.util.function.Function;

public class LNSVQDDiscretizationScheme extends MonteCarloProcessFromProcessModel {

	private final IndependentIncrements stochasticDriver;
	private final LNSVQDModel lnsvqdModel;

	/*
	 * The storage of the simulated stochastic process.
	 */
	private RandomVariable[][] discreteProcess = null;
	private transient RandomVariable[] discreteProcessWeights;

	public LNSVQDDiscretizationScheme(LNSVQDModel model, IndependentIncrements stochasticDriver) {
		super(stochasticDriver.getTimeDiscretization(), model);
		this.lnsvqdModel = model;
		this.stochasticDriver = stochasticDriver;
	}

	@Override
	public int getNumberOfPaths() {
		return stochasticDriver.getNumberOfPaths();
	}

	@Override
	public int getNumberOfFactors() {
		return stochasticDriver.getNumberOfFactors();
	}

	@Override
	public IndependentIncrements getStochasticDriver() {
		return stochasticDriver;
	}

	@Override
	public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
		return null;
	}

	@Override
	public MonteCarloProcess getCloneWithModifiedData(Map<String, Object> dataModified) {
		return null;
	}

	@Override
	public Object getCloneWithModifiedSeed(int seed) {
		return null;
	}

	/**
	 * Copied from EulerSchemeFromProcessModel
	 */
	@Override
	public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
		// Thread safe lazy initialization
		synchronized(this) {
			if (discreteProcess == null || discreteProcess.length == 0) {
				doPrecalculateProcess();
			}
		}

		if(discreteProcess[timeIndex][componentIndex] == null) {
			throw new NullPointerException("Generation of process component " + componentIndex + " at time index " + timeIndex + " failed. Likely due to out of memory");
		}

		// Return value of process
		return discreteProcess[timeIndex][componentIndex];
	}

	@Override
	public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
		return null;
	}

	@Override
	public MonteCarloProcessFromProcessModel clone() {
		return null;
	}

	/**
	 * Calculates the whole (discrete) process.
	 */
	private void doPrecalculateProcess() {
		if(discreteProcess != null && discreteProcess.length != 0) {
			return;
		}

		final int numberOfPaths = this.getNumberOfPaths();
		final int numberOfFactors = this.getNumberOfFactors();
		final int numberOfComponents = this.getNumberOfComponents();

		// Allocate Memory
		discreteProcess = new RandomVariable[getTimeDiscretization().getNumberOfTimeSteps() + 1][getNumberOfComponents()];
		discreteProcessWeights = new RandomVariable[getTimeDiscretization().getNumberOfTimeSteps() + 1];

		// Set initial Monte-Carlo weights
		discreteProcessWeights[0] = stochasticDriver.getRandomVariableForConstant(1.0 / numberOfPaths);

		// Set initial value
		final RandomVariable[] initialState = getInitialState();
		final RandomVariable[] currentState = new RandomVariable[numberOfComponents];
		for(int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++) {
			currentState[componentIndex] = initialState[componentIndex];
			discreteProcess[0][componentIndex] = applyStateSpaceTransform(0, componentIndex, currentState[componentIndex]);
		}

		// Zeta is part of the volatility discretization scheme
		Function<Double, RandomVariable> zeta = new Function<Double, RandomVariable>() {
			@Override
			public RandomVariable apply(Double aDouble) {
				return lnsvqdModel.getRandomVariableForConstant(-lnsvqdModel.getKappa1() * lnsvqdModel.getKappa2() * lnsvqdModel.getTheta() - 0.5 * lnsvqdModel.getTotalInstVar()
						+ lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() * Math.exp(aDouble) - lnsvqdModel.getKappa2()  * Math.exp(aDouble));
			}
		};

		Function<Double, RandomVariable> zeta1stDerivative = new Function<Double, RandomVariable>() {
			@Override
			public RandomVariable apply(Double aDouble) {
				return lnsvqdModel.getRandomVariableForConstant(lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() * Math.exp(aDouble) - lnsvqdModel.getKappa2()  * Math.exp(aDouble));
			}
		};

		// Evolve process
		for(int timeIndex2 = 1; timeIndex2 < getTimeDiscretization().getNumberOfTimeSteps() + 1; timeIndex2++) {
			final int timeIndex = timeIndex2;
			// Generate process from timeIndex-1 to timeIndex
			final double deltaT = getTime(timeIndex) - getTime(timeIndex - 1);

			// Fetch drift vector
			final RandomVariable[] drift;
			try {
				drift = getDrift(timeIndex - 1, discreteProcess[timeIndex - 1], null);
			} catch(final Exception e) {
				throw new RuntimeException(e + " - drift calculaton failed at time index " + timeIndex + " (time=" + getTime(timeIndex - 1) + ") . See cause of this exception for details.", e);
			}

			// Fetch brownianIncrement vector
			final RandomVariable[] brownianIncrement = stochasticDriver.getIncrement(timeIndex - 1);

			// Calculate new realization
			final ArrayList<RandomVariable> discreteProcessAtCurrentTimeIndex = new ArrayList<>(numberOfComponents);
			for(int componentIndex2 = 0; componentIndex2 < numberOfComponents; componentIndex2++) {
				final int componentIndex = componentIndex2;

				final RandomVariable driftOfComponent = drift[componentIndex];

				// Check if the component process has stopped to evolve
				if(driftOfComponent == null) {
					discreteProcessAtCurrentTimeIndex.add(componentIndex, null);
					continue;
				}

				// Apply the transform to transition into the value-space of the discretization scheme
				currentState[componentIndex] = applyStateSpaceTransform(timeIndex - 1, componentIndex, discreteProcess[timeIndex - 1][componentIndex]);
				final RandomVariable[] factorLoadings = getFactorLoading(timeIndex - 1, componentIndex, discreteProcess[timeIndex - 1]);

				/**
				 * Calculate next process value
				 */
				// Asset process
				if(componentIndex == 0 && driftOfComponent != null) {
					// Apply drift
					currentState[componentIndex] = currentState[componentIndex].sub(discreteProcess[timeIndex - 1][componentIndex].mult(-0.5 * deltaT));
					// Apply diffusion
					currentState[componentIndex] = currentState[componentIndex].addSumProduct(factorLoadings, brownianIncrement);
				}

				// Volatility process
				if(componentIndex == 1 && driftOfComponent != null) {
					double lastVolValue = discreteProcess[timeIndex - 1][componentIndex].doubleValue();
					// SOLUTION TO FORWARD ODE
					// 1. Define the function whose root is the new process-value
					Function<Double, RandomVariable> rootFunction = new Function<Double, RandomVariable>() {
						@Override
						public RandomVariable apply(Double aDouble) {
							return brownianIncrement[0].mult(-lnsvqdModel.getBeta()).sub(brownianIncrement[1].mult(lnsvqdModel.getEpsilon()))
									.sub(zeta.apply(aDouble).mult(deltaT).add(aDouble)).sub(discreteProcess[timeIndex][1]);
						}
					};

					Function<Double, RandomVariable> rootFunctionDerivative = new Function<Double, RandomVariable>() {
						@Override
						public RandomVariable apply(Double aDouble) {
							return zeta1stDerivative.apply(aDouble).mult(deltaT).sub(1);
						}
					};

					// 2. Solve for every realization
					double[] realizationsLastTimePoint = discreteProcess[timeIndex - 1][1].getRealizations();
					double[] realizationsCurrentTimePoint =  new double[getNumberOfPaths()];
					for(int j = 0; j < getNumberOfPaths(); j++) {
						double initialGuess = realizationsLastTimePoint[j]; //TODO: Make better guess
						NewtonsMethod newtonsMethod = new NewtonsMethod(initialGuess);
						// Solve
						while(!newtonsMethod.isDone()) {
							double value = rootFunction.apply(newtonsMethod.getBestPoint()).getRealizations()[j];
							double derivative = rootFunctionDerivative.apply(newtonsMethod.getBestPoint()).getRealizations()[j];
							newtonsMethod.setValueAndDerivative(value, derivative);
						}
						realizationsCurrentTimePoint[j] = newtonsMethod.getBestPoint();
					}

					currentState[componentIndex] = lnsvqdModel.getRandomVariableForArray(realizationsCurrentTimePoint);
				}

				// Transform the state space to the value space and return it.
				RandomVariable result = applyStateSpaceTransformInverse(timeIndex, componentIndex, currentState[componentIndex]);

				// The following line will add the result of the calculation to the vector discreteProcessAtCurrentTimeIndex
				discreteProcessAtCurrentTimeIndex.add(componentIndex, result);
			}

			// Fetch results and move to discreteProcess[timeIndex]
			for(int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++) {
				final RandomVariable discreteProcessAtCurrentTimeIndexAndComponent = discreteProcessAtCurrentTimeIndex.get(componentIndex);
				if(discreteProcessAtCurrentTimeIndexAndComponent != null) {
					discreteProcess[timeIndex][componentIndex] = discreteProcessAtCurrentTimeIndexAndComponent;
				} else {
					discreteProcess[timeIndex][componentIndex] = discreteProcess[timeIndex - 1][componentIndex];
				}
			}
			// Set Monte-Carlo weights
			discreteProcessWeights[timeIndex] = discreteProcessWeights[timeIndex - 1];
		}
	}

}
