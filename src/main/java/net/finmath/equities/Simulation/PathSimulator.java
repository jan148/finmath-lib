package net.finmath.equities.Simulation;

import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import net.finmath.equities.models.LNSVQDUtils;
import net.finmath.functions.NormalDistribution;
import net.finmath.modelling.descriptor.AssetModelDescriptor;
import net.finmath.modelling.descriptor.HestonModelDescriptor;
import net.finmath.modelling.descriptor.LNSVQDModelDescriptor;
import org.apache.commons.math3.random.SobolSequenceGenerator;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * This class is a path simulator for the QE-scheme for the Heston model.
 * It simulates ln(S) and V, the logarithm of the asset and the variance.
 *
 * @author Jan Berger
 */
public class PathSimulator {
	public LocalDate spotDate;
	public YieldCurve discountCurve;
	public EquityForwardStructure equityForwardStructure;
	public int numberOfPaths;
	public double[] timeGrid;
	public double[] maturities;
	public double path[][][];
	public double assetPathAtMaturities[][];

	public PathSimulator(LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure
			, int numberOfPaths, double[] timeGrid, double[] maturities) {
		this.spotDate = spotDate;
		this.discountCurve = discountCurve;
		this.equityForwardStructure = equityForwardStructure;
		this.numberOfPaths = numberOfPaths;
		this.timeGrid = timeGrid;
		this.maturities = maturities;
	}

	public void precalculatePaths(int seed, Boolean saveMemory, int startingIndex, double[] startingValue, Boolean martingaleCorrection
			, String model, String mcMethod, LNSVQDModelDescriptor lnsvqdModelDescriptor, HestonModelDescriptor hestonModelDescriptor) {
		double[] timeGridFromStartingIndex = Arrays.copyOfRange(timeGrid, startingIndex, timeGrid.length);
		double[] maturitiesFromStartingIndex = Arrays.stream(maturities).filter(x -> x >= timeGridFromStartingIndex[0]).toArray();

		BrownianBridgeNew brownianBridge = null;
		SobolSequenceGenerator sobolSequenceGenerator = null;
		if(mcMethod == "QMC") {
			ArrayList<Double> timeGridList = Arrays.stream(timeGridFromStartingIndex)
					.boxed()
					.collect(Collectors.toCollection(ArrayList::new));
			int[] prioritizedIndices = new int[maturitiesFromStartingIndex.length - 1];
			for(int j = 0; j < prioritizedIndices.length; j++) {
				double maturity = maturitiesFromStartingIndex[j];
				prioritizedIndices[j] = timeGridList.indexOf(maturity);
			}
			prioritizedIndices = prioritizedIndices.length > 0 ? prioritizedIndices : null;
			final int[][] schedulingArray = LNSVQDUtils.createSchedulingArray(timeGridFromStartingIndex.length, prioritizedIndices);

			brownianBridge = new BrownianBridgeNew(timeGridList, schedulingArray, numberOfPaths);

			int rngDim = 2 * (schedulingArray.length - 1);
			sobolSequenceGenerator = new SobolSequenceGenerator(Math.min(rngDim, 1000));
			sobolSequenceGenerator.nextVector(); // Skip 0
		}

		net.finmath.randomnumbers.MersenneTwister mersenneTwister = new net.finmath.randomnumbers.MersenneTwister(seed);

		if(assetPathAtMaturities == null) {
			assetPathAtMaturities = new double[maturities.length][numberOfPaths];
		}

		if(!saveMemory) {
			path = new double[2][timeGrid.length][numberOfPaths];
			Arrays.fill(path[0][startingIndex], Math.log(startingValue[0]));
			Arrays.fill(path[1][startingIndex], startingValue[1]);
		}

		for(int j = 0; j < numberOfPaths; j++) {
			double[][] increments = new double[2][timeGridFromStartingIndex.length - 1];
			if(mcMethod == "MC") {
				for (int r = 0; r < 2; r++) {
					for (int c = 0; c < increments[0].length; c++) {
						double deltaT = timeGridFromStartingIndex[c + 1] - timeGridFromStartingIndex[c];
						increments[r][c] = Math.sqrt(deltaT) * NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());
					}
				}
			} else {
				Random random = new Random(seed);
				double scrambleNumber = random.nextDouble();
				assert (0. < scrambleNumber && scrambleNumber < 1.) : "ScrambleNumber is out of bounds!";
				double[] vec = sobolSequenceGenerator.nextVector();
				double[] standardNormals = LNSVQDUtils.getStdNormalsFromUnifVec(vec, scrambleNumber);
				increments = brownianBridge.generateBrownianIncrementsOnePath(standardNormals, mersenneTwister);
			}

			// Fill Paths
			double asset = Math.log(startingValue[0]);
			double vol = startingValue[1];
			int currentMaturityIndex = 0;
			for(int i = 1; i < timeGridFromStartingIndex.length; i++) {
				int currentIncrementIndex = i - 1;

				double deltaT = timeGridFromStartingIndex[i] - timeGridFromStartingIndex[i - 1];
				double sqrtDeltaT = Math.sqrt(deltaT);
				assert (sqrtDeltaT > 0) : "sqrt(delta) = 0!";

				double[] incsForTimeStep = new double[]{increments[currentIncrementIndex][0], increments[currentIncrementIndex][1]};

				double[] nextVal = model == "Heston" ? hestonGetNextObs(new double[]{asset, vol}, deltaT, incsForTimeStep, hestonModelDescriptor)
						: lnsvqdGetNextObs(new double[]{asset, vol}, deltaT, incsForTimeStep, lnsvqdModelDescriptor);

				if(maturitiesFromStartingIndex[currentMaturityIndex] == timeGridFromStartingIndex[i]) {
					assetPathAtMaturities[currentMaturityIndex + (maturities.length - maturitiesFromStartingIndex.length)][j] = asset;
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

		if(martingaleCorrection) {
			for(int m = 0; m < maturities.length; m++) {
				double avg = Math.log(Arrays.stream(assetPathAtMaturities[m]).map(x -> Math.exp(x)).average().getAsDouble());
				for(int p = 0; p < numberOfPaths; p++) {
					assetPathAtMaturities[m][p] -= avg;
				}
			}
		}
	}

	private double[] hestonGetNextObs(double[] currentValue, double deltaT, double[] brownianIncrement, HestonModelDescriptor hestonModelDescriptor){
		double gamma1 = 1;
		double gamma2 = 0;
		double psiCritical = 1.5;

		double theta = hestonModelDescriptor.getTheta();
		double kappa = hestonModelDescriptor.getKappa();
		double epsilon = hestonModelDescriptor.getXi();
		double rho = hestonModelDescriptor.getRho();

		double sqrtDeltaT =  Math.sqrt(deltaT);

		// Vol
		double vol;
		double volPrev = currentValue[1]; // Needed for asset
		double m = theta + (volPrev - theta) * Math.exp(-kappa * deltaT);
		double s2 = volPrev * epsilon * epsilon * Math.exp(-kappa * deltaT) / kappa * (1 - Math.exp(-kappa * deltaT))
				+ theta * epsilon * epsilon / (2 * kappa) * (1 - Math.exp(-kappa * deltaT)) * (1 - Math.exp(-kappa * deltaT));
		double psi = s2 / (m * m);
		if(psi <= psiCritical) {
			double incVol = 1 / sqrtDeltaT * brownianIncrement[0];
			double b = Math.sqrt(2 / psi - 1 + Math.sqrt(2 / psi) * Math.sqrt(2 / psi - 1));
			double a = m / (1 + b * b);
			vol = a * (b + incVol) * (b + incVol);
		} else {
			double incVol = NormalDistribution.cumulativeDistribution(1 / sqrtDeltaT * brownianIncrement[0]) ;
			double p = (psi - 1) / (psi + 1);
			double beta = (1 - p) / m;
			double bigPsiInverse = 0 <= incVol && incVol <= p ? 0 : 1 / beta * Math.log((1 - p) / (1 - incVol));
			vol = bigPsiInverse;
		}

		// Asset
		double asset;
		double assetPrev = currentValue[0]; // Needed for asset
		double k0 = -rho * kappa * theta / epsilon * deltaT;
		double k1 = gamma1 * deltaT * (kappa * rho / epsilon - 0.5) - rho / epsilon;
		double k2 = gamma2 * deltaT * (kappa * rho / epsilon - 0.5) + rho / epsilon;
		double k3 = gamma1 * deltaT * (1 - rho * rho);
		double k4 = gamma2 * deltaT * (1 - rho * rho);
		// double A = k2 + 0.5 * k4;
		asset = assetPrev + k0 + k1 * volPrev + k2 * vol + Math.sqrt(k3 * volPrev + k4 * vol) * brownianIncrement[1];

		return new double[] {asset, vol};
	}

	private double[] lnsvqdGetNextObs(double[] currentValue, double deltaT, double[] brownianIncrement, LNSVQDModelDescriptor lnsvqdModelDescriptor){
		double kappa1 = lnsvqdModelDescriptor.getKappa1();
		double kappa2 = lnsvqdModelDescriptor.getKappa2();
		double theta = lnsvqdModelDescriptor.getTheta();
		double beta = lnsvqdModelDescriptor.getBeta();
		double epsilon = lnsvqdModelDescriptor.getEpsilon();
		double totalInstVar = lnsvqdModelDescriptor.getTotalInstVar();

		// Vol
		double vol;
		double volPrev = currentValue[1]; // Needed for asset
		double volPrevTransformed = Math.log(currentValue[1]);
		double volTransformed = volPrevTransformed + ((kappa1 * theta / volPrev - kappa1)
				+ kappa2 * (theta - volPrev) - 0.5 * totalInstVar) * deltaT
				+ beta * brownianIncrement[0]
				+ epsilon * brownianIncrement[1];
		vol = Math.exp(volTransformed);

		// Asset
		double asset;
		double assetPrev = currentValue[0];
		asset = assetPrev + volPrev * volPrev * (-0.5) * deltaT + volPrev * brownianIncrement[0];

		return new double[] {asset, vol};
	}

}
