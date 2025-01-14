package net.finmath.montecarlo.process;

import net.finmath.equities.models.LNSVQD.LNSVQDModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.rootfinder.NewtonsMethod;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.*;
import java.util.function.Function;

public class LNSVQDDiscretizationScheme extends MonteCarloProcessFromProcessModel {
	boolean newtonsMetodMultiThreaded = true;

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
	public LNSVQDModel getModel() {
		return this.lnsvqdModel;
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
			if(discreteProcess == null || discreteProcess.length == 0) {
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
		// Thread safe lazy initialization
		synchronized(this) {
			if(discreteProcessWeights == null || discreteProcessWeights.length == 0) {
				doPrecalculateProcess();
			}
		}

		// Return value of process
		return discreteProcessWeights[timeIndex];
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
			discreteProcess[0][componentIndex] = currentState[componentIndex];
		}

		// Evolve process
		for(int timeIndex2 = 1; timeIndex2 < getTimeDiscretization().getNumberOfTimeSteps() + 1; timeIndex2++) {
			// System.out.println("Currene time: " + getTime(timeIndex2));
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
					// Apply drift; Note that the drift is the drift of the log-transformed discounted process
					currentState[componentIndex] = currentState[componentIndex].add(driftOfComponent.mult(deltaT));
					// Apply diffusion; Note that the drift is the drift of the log-transformed discounted process
					currentState[componentIndex] = currentState[componentIndex].addSumProduct(factorLoadings, brownianIncrement);
				}

				// Volatility process
				if(componentIndex == 1) { /* && driftOfComponent != null*/
					// SOLUTION TO FORWARD ODE
					// 1. Define the function whose root is the new process-value
					Function<RandomVariable, RandomVariable> rootFunction = new Function<RandomVariable, RandomVariable>() {
						@Override
						public RandomVariable apply(RandomVariable randomVariable) {
							return brownianIncrement[0].mult(-lnsvqdModel.getBeta()).sub(brownianIncrement[1].mult(lnsvqdModel.getEpsilon()))
									.sub(lnsvqdModel.zeta.apply(randomVariable).mult(deltaT)).sub(currentState[componentIndex]).add(randomVariable);
							// Former: brownianIncrement[0].mult(-lnsvqdModel.getBeta()).sub(brownianIncrement[1].mult(lnsvqdModel.getEpsilon()))
							//									.sub(lnsvqdModel.zeta.apply(randomVariable).mult(deltaT)).sub(currentState[componentIndex]).add(randomVariable);
							/*return (brownianIncrement[1].mult(-Math.sqrt(lnsvqdModel.getTotalInstVar())))
									.sub(lnsvqdModel.zeta.apply(randomVariable).mult(deltaT)).sub(currentState[componentIndex]).add(randomVariable);*/
						}
					};

					Function<RandomVariable, RandomVariable> rootFunctionDerivative = new Function<RandomVariable, RandomVariable>() {
						@Override
						public RandomVariable apply(RandomVariable randomVariable) {
							return lnsvqdModel.zeta1stDerivative.apply(randomVariable).mult(deltaT).mult(-1).add(1);
						}
					};

					// 2. Solve for every realization
					double[] realizationsLastTimePoint;
					if(currentState[componentIndex] instanceof Scalar) {
						realizationsLastTimePoint = new double[getNumberOfPaths()];
						Arrays.fill(realizationsLastTimePoint, currentState[componentIndex].doubleValue());
					} else {
						realizationsLastTimePoint = currentState[componentIndex].getRealizations();
					}
					double[] realizationsCurrentTimePoint = new double[getNumberOfPaths()];
					for(int j = 0; j < getNumberOfPaths(); j++) {
						int index = j;

						/*double initialGuess = realizationsLastTimePoint[index];
						NewtonsMethod newtonsMethod = new NewtonsMethod(initialGuess);
						// Solve
						while(true) {
							double value = rootFunction.apply(lnsvqdModel.getRandomVariableForConstant(newtonsMethod.getNextPoint())).get(index); //.getRealizations()[j];
							double derivative = rootFunctionDerivative.apply(lnsvqdModel.getRandomVariableForConstant(newtonsMethod.getNextPoint())).get(index); //getRealizations()[j];
							newtonsMethod.setValueAndDerivative(value, derivative);
							if(Math.abs(value) < 1E-7) {
								break;
							}
						}
						realizationsCurrentTimePoint[index] = newtonsMethod.getBestPoint();*/
						// vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + beta*w0+volvol*w1
						realizationsCurrentTimePoint[index] = currentState[componentIndex].get(index) +
								((lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() / discreteProcess[timeIndex - 1][componentIndex].get(index) - lnsvqdModel.getKappa1())
										+ lnsvqdModel.getKappa2() * (lnsvqdModel.getTheta() - discreteProcess[timeIndex - 1][componentIndex].get(index)) - 0.5 * lnsvqdModel.getTotalInstVar()) * deltaT +
								lnsvqdModel.getBeta() * brownianIncrement[0].get(index) + lnsvqdModel.getEpsilon() * brownianIncrement[1].get(index);
						// if(realizationsCurrentTimePoint[index] < -5) {throw new RuntimeException("STOP");};
					}
					currentState[componentIndex] = lnsvqdModel.getRandomVariableForArray(realizationsCurrentTimePoint);
				}
				// Transform the state space to the value space and return it.
				currentState[componentIndex] = applyStateSpaceTransformInverse(timeIndex, componentIndex, currentState[componentIndex]);
				RandomVariable result = currentState[componentIndex];
				discreteProcessAtCurrentTimeIndex.add(componentIndex, result);
			} // End componentIndex-loop

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
