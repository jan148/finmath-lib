package net.finmath.equities.Simulation.LNSVQDPathSimulator;

import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQDModel;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.MersenneTwister;

import java.time.LocalDate;
import java.util.Arrays;

public class LNSVQDPathSimulatorMC extends LNSVQDPathSimulator{

	public LNSVQDPathSimulatorMC(LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure, int numberOfPaths, double[] timeGrid, double[] maturities, LNSVQDModel lnsvqdModel, Boolean isBackwardEuler) {
		super(spotDate, discountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturities, lnsvqdModel, isBackwardEuler);
	}

	@Override
	public void precalculatePaths(int seed, Boolean saveMemory) {
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
									- (zeta.value(Math.exp(l)) * deltaT) - volNewTransformed[pathIndex] + l)
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
					// TODO: Check replacement by zeta
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
			// TODO: Check
			// Apply martingale correction and increment maturity index
			if(maturities[currentMaturityIndex] == timeGrid[i]) {
				double avg = Math.log(Arrays.stream(assetPathAtMaturities[currentMaturityIndex]).map(x -> Math.exp(x)).average().getAsDouble());
				for(int p = 0; p < numberOfPaths; p++) {
					assetPath[p] -= avg; //+= Math.log(equityForwardStructure.getSpot() / avg);
					assetPathAtMaturities[currentMaturityIndex][p] = assetPath[p];
				}
				currentMaturityIndex++;
			}
		}
	}

}
