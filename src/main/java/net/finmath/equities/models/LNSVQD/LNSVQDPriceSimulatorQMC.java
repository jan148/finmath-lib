package net.finmath.equities.models.LNSVQD;

import net.finmath.montecarlo.*;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;

public class LNSVQDPriceSimulatorQMC extends LNSVQDEuropeanPriceSimulator {
	public LNSVQDPriceSimulatorQMC(LNSVQDModel lnsvqdModel, int numberOfPaths, double[] timeGrid, Boolean isBackwardEuler) {
		super(lnsvqdModel, numberOfPaths, timeGrid, isBackwardEuler);
	}

	public void precalculatePaths(int seed) {
		final int[][] schedulingArray = LNSVQDUtils.createSchedulingArray(timeGrid.length);

		ArrayList<Double> timeGridList = Arrays.stream(timeGrid)
				.boxed()
				.collect(Collectors.toCollection(ArrayList::new));

		/*LNSVQDUtils.printArray(timeDiscretization.getAsDoubleArray());
		System.out.println("----");
		LNSVQDUtils.printArray(timeGrid);*/
		BrownianBridgeNew brownianBridge = new BrownianBridgeNew(timeGridList, schedulingArray, numberOfPaths);

		BrentOptimizer brentOptimizer = new BrentOptimizer(1e-8, 1e-8);

		double spotPathAt0[] = new double[numberOfPaths];
		double volPathAt0[] = new double[numberOfPaths];
		Arrays.fill(spotPathAt0, Math.log(lnsvqdModel.getSpot0()));
		Arrays.fill(volPathAt0, lnsvqdModel.getSigma0());

		path[0][0] = spotPathAt0;
		path[1][0] = volPathAt0;
		double[] volNewTransformed = new double[numberOfPaths];
		Arrays.fill(volNewTransformed, Math.log(lnsvqdModel.getSigma0()));

		for(int i = 1; i < timeGrid.length; i++) {
			double deltaT = timeGrid[i] - timeGrid[i - 1];
			double sqrtDeltaT = Math.sqrt(deltaT);
			assert(sqrtDeltaT > 0) : "sqrt(delta) = 0!";
			double[][] brownianIncrements = new double[numberOfPaths][2];
			// Fill Paths
			for(int j = 0; j < numberOfPaths; j++) {
				int pathIndex = j;

				brownianIncrements[j][0] = brownianBridge.getBrownianIncrementArr(i - 1, 0, seed)[j];
				brownianIncrements[j][1] = brownianBridge.getBrownianIncrementArr(i - 1, 1, seed)[j];

				if(isBackwardEuler) {
					UnivariateObjectiveFunction rootFunction = new UnivariateObjectiveFunction(
							l -> Math.abs(-brownianIncrements[pathIndex][0] * lnsvqdModel.getBeta() - (brownianIncrements[pathIndex][1] * lnsvqdModel.getEpsilon())
									- (zeta.value(l) * deltaT) - volNewTransformed[pathIndex] + l)
					);
					UnivariatePointValuePair result = brentOptimizer.optimize(
							rootFunction,
							GoalType.MINIMIZE,
							new MaxEval(100),
							new org.apache.commons.math3.optim.univariate.SearchInterval(-100, 100)
					);
					volNewTransformed[j] = result.getPoint();
					if(Math.abs(result.getValue()) > 1e-4) {throw new ArithmeticException("The point doesn't result in a root.");}
				} else {
					// from Sepp's implementation: vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + vartheta*w1_
					//				sigma0 = np.exp(vol_var)
					/*sigma0_2dt = vol_backbone_eta2 * sigma0 * sigma0 * dt
					x0 = x0 + alpha * 0.5 * sigma0_2dt + vol_backbone_eta * sigma0 * w0
					vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + beta*w0+volvol*w1
					sigma0 = np.exp(vol_var)*/
					volNewTransformed[j] = volNewTransformed[j] + ((lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() / path[1][i - 1][j] - lnsvqdModel.getKappa1())
							+ lnsvqdModel.getKappa2() * (lnsvqdModel.getTheta() - path[1][i - 1][j]) - 0.5 * lnsvqdModel.getTotalInstVar()) * deltaT
							+ lnsvqdModel.getBeta() * brownianIncrements[j][0] + lnsvqdModel.getEpsilon() * brownianIncrements[j][1];
				}
				// vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + vartheta*w1_
				path[1][i][j] = Math.exp(volNewTransformed[j]);

				// Asset path
				path[0][i][j] = path[0][i - 1][j] + path[1][i - 1][j] * path[1][i - 1][j] * (-0.5) * deltaT + path[1][i - 1][j] * brownianIncrements[j][0];
			}
		}
	}

}
