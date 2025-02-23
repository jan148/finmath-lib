package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.models.Black76Model;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.MersenneTwister;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class LNSVQDCallPriceSimulator {
	LNSVQDModel lnsvqdModel;
	int numberOfPaths;
	double[] timeGrid;
	double path[][][];
	boolean isBackwardEuler = false;
	UnivariateFunction zeta = x -> Math.exp(-x) * lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() - Math.exp(x) * lnsvqdModel.getKappa2()
			- lnsvqdModel.getKappa1() + lnsvqdModel.getKappa2() * lnsvqdModel.getTheta() - 0.5 * lnsvqdModel.getTotalInstVar();

	public LNSVQDCallPriceSimulator(LNSVQDModel lnsvqdModel, int numberOfPaths, double[] timeGrid, Boolean isBackwardEuler) {
		this.lnsvqdModel = lnsvqdModel;
		this.numberOfPaths = numberOfPaths;
		this.timeGrid = timeGrid;
		// Component, time, paths
		this.path = new double[2][timeGrid.length][numberOfPaths];
		this.isBackwardEuler = isBackwardEuler;
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
			double discountFactor = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(timeGrid[i - 1]);
			double discountFactorCurrent = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(timeGrid[i]);
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

				if(isBackwardEuler) {
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
					// from Sepp's implementation: vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + vartheta*w1_
					//				sigma0 = np.exp(vol_var)
					/*sigma0_2dt = vol_backbone_eta2 * sigma0 * sigma0 * dt
					x0 = x0 + alpha * 0.5 * sigma0_2dt + vol_backbone_eta * sigma0 * w0
					vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + beta*w0+volvol*w1
					sigma0 = np.exp(vol_var)*/
					volNewTransformed = volPrevTransformed + ((lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() / volPrev - lnsvqdModel.getKappa1())
							+ lnsvqdModel.getKappa2() * (lnsvqdModel.getTheta() - volPrev) - 0.5 * lnsvqdModel.getTotalInstVar()) * deltaT
							+ lnsvqdModel.getBeta() * brownianIncrements[j][0] + lnsvqdModel.getEpsilon() * brownianIncrements[j][1];
				}
				// vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + vartheta*w1_
				path[1][i][j] = Math.exp(volNewTransformed);

				// Vol path
				double assetPrev = path[0][i - 1][j];
				double assetTransformed = Math.log(assetPrev * discountFactor);
				assetTransformed += Math.pow(volPrev, 2) * (-0.5) * deltaT + volPrev * brownianIncrements[j][0];
				path[0][i][j] = Math.exp(assetTransformed) / discountFactorCurrent;
			}
		}
	}

	public void precalculatePathsNew(int seed) {
		MersenneTwister mersenneTwister = new MersenneTwister(seed);
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

				brownianIncrements[j][0] = mersenneTwister.nextGaussian() * sqrtDeltaT;
				brownianIncrements[j][1] = mersenneTwister.nextGaussian() * sqrtDeltaT;

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

	public double[] getPayoffsAtMaturity(double strike, double maturity, int callPutSign) throws Exception {
		List<Double> list = Arrays.stream(timeGrid).boxed().collect(Collectors.toList());
		int matIndex = list.indexOf(maturity);
		if(matIndex == -1) {throw new Exception("Maturity not found!");}
		double forward = lnsvqdModel.equityForwardStructure.getForward(maturity) * lnsvqdModel.getSpot0();
		// In Sepp: spots_t = forward*np.exp(x0)
		//				correnction = np.nanmean(spots_t) - forward
		//				spots_t = spots_t - correnction
		double[] actualAssets = Arrays.stream(path[0][matIndex]).map(x -> Math.exp(x) * forward).toArray();
		for(double a : actualAssets) {assert(!Double.isNaN(a)) : "Nan encountered";}
		double correction = Arrays.stream(actualAssets).average().getAsDouble() - forward;
		actualAssets = Arrays.stream(actualAssets).map(x -> x - correction).toArray();
		double[] payoffsAtMaturity = Arrays.stream(actualAssets)
				.map(x -> Math.max(callPutSign * (x - strike), 0)).toArray();
		return payoffsAtMaturity;
	}

	public double getCallPrice(double strike, double maturity, int callPutSign) throws Exception {
		double discountFactor = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
		double[] payoffsAtMaturity = getPayoffsAtMaturity(strike, maturity, callPutSign);
		for(double payoff : payoffsAtMaturity) {assert(!Double.isNaN(payoff)) : "Nan encountered";}
		double expectationAtMaturity = Arrays.stream(payoffsAtMaturity).average().getAsDouble();
		double price = expectationAtMaturity * discountFactor;
		return price;
	}

	public double getImpliedVol(double strike, double maturity) throws Exception {
		double discountFactor = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
		double forward = lnsvqdModel.equityForwardStructure.getForward(maturity) * lnsvqdModel.getSpot0();

		int callPutSign;
		if(strike > forward) {
			callPutSign = 1;
		} else {
			callPutSign = -1;
		}
		double price = getCallPrice(strike, maturity, callPutSign);

		double impliedVol;
		if(strike > forward) {
			impliedVol = Black76Model.optionImpliedVolatility(forward, strike, maturity, price / discountFactor, true);
		} else {
			impliedVol = Black76Model.optionImpliedVolatility(forward, strike, maturity, price / discountFactor, false);
		}
		return impliedVol;
	}
}
