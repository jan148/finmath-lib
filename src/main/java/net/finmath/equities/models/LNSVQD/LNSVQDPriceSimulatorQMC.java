package net.finmath.equities.models.LNSVQD;

import net.finmath.montecarlo.BrownianBridge;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromRandomNumberGenerator;
import net.finmath.randomnumbers.SobolSequence;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.SobolSequenceGenerator;

import java.util.Arrays;

public class LNSVQDPriceSimulatorQMC extends LNSVQDCallPriceSimulator{

	public LNSVQDPriceSimulatorQMC(LNSVQDModel lnsvqdModel, int numberOfPaths, double[] timeGrid) {
		super(lnsvqdModel, numberOfPaths, timeGrid);
	}

	public void precalculatePaths(int seed) {
		TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(timeGrid);
		// Dim = 2 * (timeGrid - 1)
		SobolSequence sobolSequenceGenerator = new SobolSequence(2 * (timeGrid.length - 1));
		sobolSequenceGenerator.generator.skipTo(Math.abs(seed) * numberOfPaths);

		BrownianMotion generator = new BrownianMotionFromRandomNumberGenerator(timeDiscretization, 2, numberOfPaths, sobolSequenceGenerator);

		RandomVariable[] start = new RandomVariable[]{generator.getRandomVariableForConstant(0), generator.getRandomVariableForConstant(0)};
		TimeDiscretization timeDiscretizationForEnd = new TimeDiscretizationFromArray(timeGrid[0], timeGrid[timeGrid.length - 1]);
		// Dim = 2 * (timeGrid - 1)
		SobolSequence sobolSequenceGeneratorForEnd = new SobolSequence(2);
		sobolSequenceGeneratorForEnd.generator.skipTo(Math.abs(seed) * numberOfPaths);
		BrownianMotion generatorForEnd = new BrownianMotionFromRandomNumberGenerator(timeDiscretizationForEnd, 2, numberOfPaths, sobolSequenceGeneratorForEnd);
		RandomVariable[] end = new RandomVariable[]{generatorForEnd.getBrownianIncrement(0, 0), generatorForEnd.getBrownianIncrement(0, 1)};
		BrownianBridge brownianBridge = new BrownianBridge(generator, start, end);

		BrentOptimizer brentOptimizer = new BrentOptimizer(1e-8, 1e-8);

		double spotPathAt0[] = new double[numberOfPaths];
		double volPathAt0[] = new double[numberOfPaths];
		Arrays.fill(spotPathAt0, lnsvqdModel.getSpot0());
		Arrays.fill(volPathAt0, lnsvqdModel.getSigma0());

		path[0][0] = spotPathAt0;
		path[1][0] = volPathAt0;

		for(int i = 1; i < timeGrid.length; i++) {
			double deltaT = timeGrid[i] - timeGrid[i - 1];
			// double sqrtDeltaT = Math.sqrt(deltaT);
			double discountFactor = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(timeGrid[i - 1]);
			double discountFactorCurrent = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(timeGrid[i]);
			double[][] brownianIncrements = new double[numberOfPaths][2];
			// Fill Paths
			for(int j = 0; j < numberOfPaths; j++) {
				int pathIndex = j;

				brownianIncrements[j][0] = brownianBridge.getBrownianIncrement(i - 1, 0).get(j);
				brownianIncrements[j][1] =brownianBridge.getBrownianIncrement(i - 1, 1).get(j);

				// Vol path
				double volPrev = path[1][i - 1][j];
				double volPrevTransformed = Math.log(volPrev);
				double volNewTransformed;

				if(isForwardEuler) {
					UnivariateObjectiveFunction rootFunction = new UnivariateObjectiveFunction(
							l -> Math.abs(-brownianIncrements[pathIndex][0] * lnsvqdModel.getBeta() - (brownianIncrements[pathIndex][1] * lnsvqdModel.getEpsilon())
									- (zeta.value(l) * deltaT) - volPrevTransformed + l)
					);
					UnivariatePointValuePair result = brentOptimizer.optimize(
							rootFunction,
							GoalType.MINIMIZE,
							new MaxEval(100),
							new org.apache.commons.math3.optim.univariate.SearchInterval(-1000, 1000)
					);
					volNewTransformed = result.getPoint();
					if(Math.abs(result.getValue()) > 1e-4) {throw new ArithmeticException("The point doesn't result in a root.");}
				} else {
					volNewTransformed = volPrevTransformed + ((lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() / volPrev - lnsvqdModel.getKappa1())
							+ lnsvqdModel.getKappa2() * (lnsvqdModel.getTheta() - volPrev) - 0.5 * lnsvqdModel.getTotalInstVar()) * deltaT +
							lnsvqdModel.getBeta() * brownianIncrements[j][0] + lnsvqdModel.getEpsilon() * brownianIncrements[j][1];
				}
				path[1][i][j] = Math.exp(volNewTransformed);

				// Vol path
				double assetPrev = path[0][i - 1][j];
				double assetTransformed = Math.log(assetPrev * discountFactor);
				assetTransformed += Math.pow(volPrev, 2) * (-0.5) * deltaT + volPrev * brownianIncrements[j][0];
				path[0][i][j] = Math.exp(assetTransformed) / discountFactorCurrent;
			}
		}
	}
}
