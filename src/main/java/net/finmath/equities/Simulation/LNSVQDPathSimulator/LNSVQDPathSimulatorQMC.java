package net.finmath.equities.Simulation.LNSVQDPathSimulator;

import net.finmath.equities.Simulation.BrownianBridgeNew;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.equities.models.LNSVQDModel;
import net.finmath.randomnumbers.MersenneTwister;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.SobolSequenceGenerator;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.Collectors;

public class LNSVQDPathSimulatorQMC extends LNSVQDPathSimulator {

	public LNSVQDPathSimulatorQMC(LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure, int numberOfPaths, double[] timeGrid, double[] maturities, LNSVQDModel lnsvqdModel, Boolean isBackwardEuler) {
		super(spotDate, discountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturities, lnsvqdModel, isBackwardEuler);
	}

	@Override
	public void precalculatePaths(int seed, Boolean saveMemory, int startingIndex, double[] startingValue, Boolean martingaleCorrection) {
		ArrayList<Double> timeGridList = Arrays.stream(timeGrid)
				.boxed()
				.collect(Collectors.toCollection(ArrayList::new));
		int[] prioritizedIndices = new int[maturities.length - 1];
		for(int j = 0; j < prioritizedIndices.length; j++) {
			double maturity = maturities[j];
			prioritizedIndices[j] = timeGridList.indexOf(maturity);
		}
		prioritizedIndices = prioritizedIndices.length > 0 ? prioritizedIndices : null;
		final int[][] schedulingArray = LNSVQDUtils.createSchedulingArray(timeGrid.length, prioritizedIndices);

		BrownianBridgeNew brownianBridge = new BrownianBridgeNew(timeGridList, schedulingArray, numberOfPaths);

		/**
		 * Init Sobol and Mersenne
		 */
		int rngDim = 2 * (schedulingArray.length - 1);
		SobolSequenceGenerator sobolSequenceGenerator = new SobolSequenceGenerator(Math.min(rngDim, 1000));
		sobolSequenceGenerator.nextVector(); // Skip 0

		Random random = new Random(seed);
		double scrambleNumber = random.nextDouble();
		assert (0. < scrambleNumber && scrambleNumber < 1.) : "ScrambleNumber is out of bounds!";

		MersenneTwister mersenneTwister = new MersenneTwister(seed);
		/**
		 * Optimizer
		 */
		BrentOptimizer brentOptimizer = new BrentOptimizer(1e-8, 1e-8);

		if(assetPathAtMaturities == null) {
			assetPathAtMaturities = new double[maturities.length][numberOfPaths];
		}

		if(!saveMemory) {
			path = new double[2][timeGrid.length][numberOfPaths];
			Arrays.fill(path[0][0], Math.log(startingValue[0]));
			Arrays.fill(path[1][0], startingValue[1]);
		}

		for(int j = 0; j < numberOfPaths; j++) {
			double[] vec = sobolSequenceGenerator.nextVector();
			double[] standardNormals = LNSVQDUtils.getStdNormalsFromUnifVec(vec, scrambleNumber);
			double[][] brownianIncrements = brownianBridge.generateBrownianIncrementsOnePath(standardNormals, mersenneTwister); // new double[numberOfPaths][2];
			// Fill Paths
			double asset = Math.log(startingValue[0]);
			double vol = startingValue[1];
			double volTransformed = Math.log(startingValue[1]);
			int currentMaturityIndex = Arrays.stream(maturities).filter(x -> x < timeGrid[startingIndex]).toArray().length;
			for(int i = startingIndex; i < timeGrid.length; i++) {
				int currentIncrementIndex = i - 1;

				double deltaT = timeGrid[i] - timeGrid[i - 1];
				if(isBackwardEuler) {
					double copyVolTransformed = volTransformed; // Need to copy bc. of static context
					UnivariateObjectiveFunction rootFunction = new UnivariateObjectiveFunction(
							l -> -brownianIncrements[currentIncrementIndex][1] * lnsvqdModel.getBeta() - (brownianIncrements[currentIncrementIndex][0] * lnsvqdModel.getEpsilon())
									- (zeta.value(Math.exp(l)) * deltaT) - copyVolTransformed + l
					);
					UnivariatePointValuePair result = brentOptimizer.optimize(
							rootFunction,
							GoalType.MINIMIZE,
							new MaxEval(100),
							new org.apache.commons.math3.optim.univariate.SearchInterval(-100, 100)
					);
					volTransformed = result.getPoint();
					if(Math.abs(result.getValue()) > 1e-4) {
						throw new ArithmeticException("The point doesn't result in a root.");
					}
				} else {
					// TODO: Check replacement by zeta
					volTransformed = volTransformed + ((lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() / vol - lnsvqdModel.getKappa1())
							+ lnsvqdModel.getKappa2() * (lnsvqdModel.getTheta() - vol) - 0.5 * lnsvqdModel.getTotalInstVar()) * deltaT
							+ lnsvqdModel.getBeta() * brownianIncrements[currentIncrementIndex][1] + lnsvqdModel.getEpsilon() * brownianIncrements[currentIncrementIndex][0];
				}
				asset = asset + vol * vol * (-0.5) * deltaT + vol * brownianIncrements[currentIncrementIndex][1];
				vol = Math.exp(volTransformed);

				// System.out.println(maturities[currentMaturityIndex] + "\t" + timeGrid[i]);
				if(maturities[currentMaturityIndex] == timeGrid[i]) {
					assetPathAtMaturities[currentMaturityIndex][j] = asset;
					currentMaturityIndex = currentMaturityIndex + 1;
				}
				if(!saveMemory) {
					// Vol path
					path[1][i][j] = vol;
					// Asset path
					path[0][i][j] = asset;
				}
			}
		}
		// Apply martingale correction; Can only apply it after complete rollout, might be problematic
		if(martingaleCorrection) {
			for(int m = 0; m < maturities.length; m++) {
				double avg = Math.log(Arrays.stream(assetPathAtMaturities[m]).map(x -> Math.exp(x)).average().getAsDouble());
				for(int p = 0; p < numberOfPaths; p++) {
					assetPathAtMaturities[m][p] -= avg;
				}
			}
		}
	}

}
