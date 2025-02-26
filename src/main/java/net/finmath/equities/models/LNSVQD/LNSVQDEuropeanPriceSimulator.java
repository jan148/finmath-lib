package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.models.Black76Model;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class LNSVQDEuropeanPriceSimulator {
	LNSVQDModel lnsvqdModel;
	int numberOfPaths;
	double[] timeGrid;
	double[] maturities;
	double path[][][];
	double assetPathAtMaturities[][];
	boolean isBackwardEuler = false;
	UnivariateFunction zeta = x -> Math.exp(-x) * lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() - Math.exp(x) * lnsvqdModel.getKappa2()
			- lnsvqdModel.getKappa1() + lnsvqdModel.getKappa2() * lnsvqdModel.getTheta() - 0.5 * lnsvqdModel.getTotalInstVar();

	public LNSVQDEuropeanPriceSimulator(LNSVQDModel lnsvqdModel, int numberOfPaths, double[] timeGrid, double[] maturities, Boolean isBackwardEuler) {
		this.lnsvqdModel = lnsvqdModel;
		this.numberOfPaths = numberOfPaths;
		this.timeGrid = timeGrid;
		this.maturities = maturities;
		// Component, time, paths
		this.isBackwardEuler = isBackwardEuler;
	}

	public void precalculatePaths(int seed, Boolean saveMemory) throws Exception {
		MersenneTwister mersenneTwister = new MersenneTwister(seed);
		BrentOptimizer brentOptimizer = new BrentOptimizer(1e-8, 1e-8);
		int currentMaturityIndex = 0;

		assetPathAtMaturities = new double[maturities.length][numberOfPaths];
		if(!saveMemory) {
			path = new double[2][timeGrid.length][numberOfPaths];
		}

		double assetPath[] = new double[numberOfPaths];
		double volPath[] = new double[numberOfPaths];
		Arrays.fill(assetPath, Math.log(lnsvqdModel.getSpot0()));
		Arrays.fill(volPath, lnsvqdModel.getSigma0());

		if(!saveMemory) {
			path[0][0] = assetPath;
			path[1][0] = volPath;
		}

		double[] volNewTransformed = new double[numberOfPaths];
		Arrays.fill(volNewTransformed, Math.log(lnsvqdModel.getSigma0()));

		for(int i = 1; i < timeGrid.length; i++) {
			double deltaT = timeGrid[i] - timeGrid[i - 1];
			double sqrtDeltaT = Math.sqrt(deltaT);
			assert (sqrtDeltaT > 0) : "sqrt(delta) = 0!";
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
					if(Math.abs(result.getValue()) > 1e-4) {
						throw new ArithmeticException("The point doesn't result in a root.");
					}
				} else {
					// from Sepp's implementation: vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + vartheta*w1_
					//				sigma0 = np.exp(vol_var)
					/*sigma0_2dt = vol_backbone_eta2 * sigma0 * sigma0 * dt
					x0 = x0 + alpha * 0.5 * sigma0_2dt + vol_backbone_eta * sigma0 * w0
					vol_var = vol_var + ((kappa1 * theta / sigma0 - kappa1) + kappa2*(theta-sigma0) + adj*sigma0 - 0.5*vartheta2) * dt + beta*w0+volvol*w1
					sigma0 = np.exp(vol_var)*/
					volNewTransformed[j] = volNewTransformed[j] + ((lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() / volPath[j] - lnsvqdModel.getKappa1())
							+ lnsvqdModel.getKappa2() * (lnsvqdModel.getTheta() - volPath[j]) - 0.5 * lnsvqdModel.getTotalInstVar()) * deltaT
							+ lnsvqdModel.getBeta() * brownianIncrements[j][0] + lnsvqdModel.getEpsilon() * brownianIncrements[j][1];
				}
				assetPath[j] = assetPath[j] + volPath[j] * volPath[j] * (-0.5) * deltaT + volPath[j] * brownianIncrements[j][0];
				volPath[j] = Math.exp(volNewTransformed[j]);

				if(maturities[currentMaturityIndex] == timeGrid[i]) {
					assetPathAtMaturities[currentMaturityIndex][j] = assetPath[j];
				}
				if(!saveMemory) {
					// Vol path
					path[1][i][j] = volPath[j];
					// Asset path
					path[0][i][j] = assetPath[j];
				}
			}
			if(maturities[currentMaturityIndex] == timeGrid[i]) {currentMaturityIndex++;}
		}
	}

	public double[] getPayoffsAtMaturity(double strike, double maturity, int callPutSign) throws Exception {
		List<Double> list = Arrays.stream(maturities).boxed().collect(Collectors.toList());
		int matIndex = list.indexOf(maturity);
		if(matIndex == -1) {
			throw new Exception("Maturity not found!");
		}
		double forward = lnsvqdModel.equityForwardStructure.getForward(maturity) / lnsvqdModel.spot0;
		double[] actualAssets = Arrays.stream(assetPathAtMaturities[matIndex]).map(x -> Math.exp(x) * forward).toArray();
		double[] payoffsAtMaturity = Arrays.stream(actualAssets)
				.map(x -> Math.max(callPutSign * (x - strike), 0)).toArray();
		return payoffsAtMaturity;
	}

	public double getEuropeanPrice(double strike, double maturity, int callPutSign) throws Exception {
		double discountFactor = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
		double[] payoffsAtMaturity = getPayoffsAtMaturity(strike, maturity, callPutSign);
		// for(double payoff : payoffsAtMaturity) {assert(!Double.isNaN(payoff)) : "Nan encountered";}
		double expectationAtMaturity = Arrays.stream(payoffsAtMaturity).average().getAsDouble();
		double price = expectationAtMaturity * discountFactor;
		return price;
	}

	public double getEuropeanPriceAuto(double strike, double maturity) throws Exception {
		double forward = lnsvqdModel.equityForwardStructure.getForward(maturity);
		int callPutSign;
		if(strike > forward) {
			callPutSign = 1;
		} else {
			callPutSign = -1;
		}
		return getEuropeanPrice(strike, maturity, callPutSign);
	}

	public double[] getPriceStdErrorAndBounds(double strike, double maturity) throws Exception {
		double[] result = new double[4];
		StandardDeviation standardDeviation = new StandardDeviation(true);
		double discountFactor = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
		double forward = lnsvqdModel.equityForwardStructure.getForward(maturity);

		int callPutSign;
		if(strike > forward) {
			callPutSign = 1;
		} else {
			callPutSign = -1;
		}

		double[] payoffsAtMaturity = getPayoffsAtMaturity(strike, maturity, callPutSign);
		for(double payoff : payoffsAtMaturity) {
			assert (!Double.isNaN(payoff)) : "Nan encountered";
		}
		double expectationAtMaturity = Arrays.stream(payoffsAtMaturity).average().getAsDouble();
		double price = expectationAtMaturity * discountFactor;
		result[0] = price;
		double stdError = standardDeviation.evaluate(payoffsAtMaturity) * discountFactor;
		result[1] = stdError;
		double[] boundsUndiscounted = LNSVQDUtils.getConfidenceInterval(payoffsAtMaturity, 0.05);
		result[2] = boundsUndiscounted[0] * discountFactor;
		result[3] = boundsUndiscounted[1] * discountFactor;
		return result;
	}

	public double getImpliedVol(double strike, double maturity) throws Exception {
		double discountFactor = lnsvqdModel.equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
		double forward = lnsvqdModel.equityForwardStructure.getForward(maturity);

		int callPutSign;
		if(strike > forward) {
			callPutSign = 1;
		} else {
			callPutSign = -1;
		}
		double price = getEuropeanPrice(strike, maturity, callPutSign);

		double impliedVol;
		if(strike > forward) {
			impliedVol = Black76Model.optionImpliedVolatility(forward, strike, maturity, price / discountFactor, true);
		} else {
			impliedVol = Black76Model.optionImpliedVolatility(forward, strike, maturity, price / discountFactor, false);
		}
		return impliedVol;
	}
}
