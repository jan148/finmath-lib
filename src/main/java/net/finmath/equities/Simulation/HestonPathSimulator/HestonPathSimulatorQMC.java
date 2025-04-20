package net.finmath.equities.Simulation.HestonPathSimulator;

import net.finmath.equities.Simulation.BrownianBridgeNew;
import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.functions.NormalDistribution;
import org.apache.commons.math3.random.SobolSequenceGenerator;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Implementation of QE-scheme for Heston; Simulation of ln X and V
 */
public class HestonPathSimulatorQMC extends HestonPathSimulator {
	double gamma1;
	double gamma2;
	double psiCritical = 1.5; // Default in paper

	public HestonPathSimulatorQMC(LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure, int numberOfPaths, double[] timeGrid, double[] maturities, double sigma0, double kappa, double theta, double epsilon, double rho, double gamma1, double gamma2) {
		super(spotDate, discountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturities, sigma0, kappa, theta, epsilon, rho);
		this.gamma1 = gamma1;
		this.gamma2 = gamma2;
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

		net.finmath.randomnumbers.MersenneTwister mersenneTwister = new net.finmath.randomnumbers.MersenneTwister(seed);

		if(assetPathAtMaturities == null) {assetPathAtMaturities = new double[maturities.length][numberOfPaths];}

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
			int currentMaturityIndex = Arrays.stream(maturities).filter(x -> x < timeGrid[startingIndex]).toArray().length;
			for(int i = startingIndex; i < timeGrid.length; i++) {
				int currentIncrementIndex = i - 1;

				double deltaT = timeGrid[i] - timeGrid[i - 1];
				double sqrtDeltaT = Math.sqrt(deltaT);
				assert (sqrtDeltaT > 0) : "sqrt(delta) = 0!";

				// Pre-chaching
				// Asset
				double k0 = -rho * kappa * theta / epsilon * deltaT;
				double k1 = gamma1 * deltaT * (kappa * rho / epsilon - 0.5) - rho / epsilon;
				double k2 = gamma2 * deltaT * (kappa * rho / epsilon - 0.5) + rho / epsilon;
				double k3 = gamma1 * deltaT * (1 - rho * rho);
				double k4 = gamma2 * deltaT * (1 - rho * rho);
				// double A = k2 + 0.5 * k4;

				// Fill Paths
				double incAsset = 1 / sqrtDeltaT * brownianIncrements[currentIncrementIndex][1];

				// Vol
				double incVol;
				double volPrev = vol; // Needed for asset
				double m = theta + (volPrev - theta) * Math.exp(-kappa * deltaT);
				double s2 = volPrev * epsilon * epsilon * Math.exp(-kappa * deltaT) / kappa * (1 - Math.exp(-kappa * deltaT))
						+ theta * epsilon * epsilon / (2 * kappa) * (1 - Math.exp(-kappa * deltaT)) * (1 - Math.exp(-kappa * deltaT));
				double psi = s2 / (m * m);
				if(psi <= psiCritical) {
					incVol = 1 / sqrtDeltaT * brownianIncrements[currentIncrementIndex][0];
					double b = Math.sqrt(2 / psi - 1 + Math.sqrt(2 / psi) * Math.sqrt(2 / psi - 1));
					double a = m / (1 + b * b);
					vol = a * (b + incVol) * (b + incVol);
				} else {
					incVol = NormalDistribution.cumulativeDistribution(1 / sqrtDeltaT * brownianIncrements[currentIncrementIndex][0]) ;
					double p = (psi - 1) / (psi + 1);
					double beta = (1 - p) / m;
					// System.out.println((1 - p) / beta + "\t" + m);
					// System.out.println((1 - p * p) / (beta * beta) + "\t" + s2);
					double bigPsiInverse = 0 <= incVol && incVol <= p ? 0 : 1 / beta * Math.log((1 - p) / (1 - incVol));
					vol = bigPsiInverse;
				}

				// Asset
				asset = asset + k0 + k1 * volPrev + k2 * vol + Math.sqrt(k3 * volPrev + k4 * vol) * incAsset;

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
		// TODO: Check
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
