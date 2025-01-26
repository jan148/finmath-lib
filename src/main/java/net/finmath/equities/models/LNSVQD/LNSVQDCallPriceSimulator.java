package net.finmath.equities.models.LNSVQD;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.MersenneTwister;

import java.util.Arrays;

public class LNSVQDCallPriceSimulator {
	LNSVQDModel lnsvqdModel;
	int numberOfPaths;
	double[] timeGrid;
	double path[][][];
	boolean isForwardEuler = false;
	UnivariateFunction zeta = x -> Math.exp(-x) * lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() - Math.exp(x) * lnsvqdModel.getKappa2()
			- lnsvqdModel.getKappa1() + lnsvqdModel.getKappa2() * lnsvqdModel.getTheta() - 0.5 * lnsvqdModel.getTotalInstVar();

	public LNSVQDCallPriceSimulator(LNSVQDModel lnsvqdModel, int numberOfPaths, double[] timeGrid) {
		this.lnsvqdModel = lnsvqdModel;
		this.numberOfPaths = numberOfPaths;
		this.timeGrid = timeGrid;
		// Component, time, paths
		this.path = new double[2][timeGrid.length][numberOfPaths];
	}

	public double[][][] getTransformedPath() {
		if(path == null) {
			throw new IllegalStateException("Path hasn't been precalculated yet!");
		} else {
			double[][][] transformedPath = this.path.clone();
			for(int i = 1; i < timeGrid.length; i++) {
				double deltaT = timeGrid[i] - timeGrid[i - 1];
				double sqrtDeltaT = Math.sqrt(deltaT);
				double discountFactor = Math.exp(-lnsvqdModel.getRiskFreeRate(timeGrid[i]) * timeGrid[i]);
				double[][] brownianIncrements = new double[numberOfPaths][2];
				// Fill Paths
				for(int j = 0; j < numberOfPaths; j++) {
					transformedPath[0][i][j] = Math.log(transformedPath[0][i][j] * discountFactor);
					transformedPath[0][i][j] = transformedPath[1][i][j] - lnsvqdModel.theta;
				}
			}
			return transformedPath;
		}
	}


	public void precalculatePaths(int seed) {
		MersenneTwister mersenneTwister = new MersenneTwister(seed);
		BrentOptimizer brentOptimizer = new BrentOptimizer(1e-8, 1e-8);

		double spotPathAt0[] = new double[numberOfPaths];
		double volPathAt0[] = new double[numberOfPaths];
		Arrays.fill(spotPathAt0, lnsvqdModel.getSpot0());
		Arrays.fill(volPathAt0, lnsvqdModel.getSigma0());

		path[0][0] = spotPathAt0;
		path[1][0] = volPathAt0;

		for(int i = 1; i < timeGrid.length; i++) {
			double deltaT = timeGrid[i] - timeGrid[i - 1];
			double sqrtDeltaT = Math.sqrt(deltaT);
			double discountFactor = Math.exp(-lnsvqdModel.getRiskFreeRate(timeGrid[i]) * deltaT);
			double[][] brownianIncrements = new double[numberOfPaths][2];
			// Fill Paths
			for(int j = 0; j < numberOfPaths; j++) {
				int pathIndex = j;

				brownianIncrements[j][0] = mersenneTwister.nextGaussian() * sqrtDeltaT;
				brownianIncrements[j][1] = mersenneTwister.nextGaussian() * sqrtDeltaT;

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
				path[0][i][j] = Math.exp(assetTransformed) / discountFactor;
			}
		}
	}

	public double getCallPrice(double strike) {
		double discountFactor = Math.exp(-lnsvqdModel.getRiskFreeRate(0) * timeGrid[timeGrid.length - 1]); // Todo: Change (not index 0)!
		double expectationAtMaturity = Arrays.stream(path[0][timeGrid.length - 1])
				.map(x -> Math.max(x - strike, 0)).average().getAsDouble();
		double price = expectationAtMaturity * discountFactor;
		return price;
	}

}
