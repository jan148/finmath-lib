package net.finmath.equities.Simulation.HestonPathSimulator;

import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;
import org.apache.commons.math3.random.MersenneTwister;

import java.time.LocalDate;
import java.util.Arrays;

/**
 * Implementation of QE-scheme for Heston; Simulation of ln X and V
 */
public class HestonPathSimulatorMC extends HestonPathSimulator {
	double gamma1;
	double gamma2;
	double psiCritical = 1.5; // Default in paper

	public HestonPathSimulatorMC(LocalDate spotDate, YieldCurve discountCurve, EquityForwardStructure equityForwardStructure, int numberOfPaths, double[] timeGrid, double[] maturities, double sigma0, double kappa, double theta, double epsilon, double rho, double gamma1, double gamma2) {
		super(spotDate, discountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturities, sigma0, kappa, theta, epsilon, rho);
		this.gamma1 = gamma1;
		this.gamma2 = gamma2;
	}

	@Override
	public void precalculatePaths(int seed, Boolean saveMemory) {
		MersenneTwister mersenneTwister = new MersenneTwister(seed);
		int currentMaturityIndex = 0;

		assetPathAtMaturities = new double[maturities.length][numberOfPaths];
		if(!saveMemory) {
			path = new double[2][timeGrid.length][numberOfPaths];
		}

		double assetPath[] = new double[numberOfPaths];
		double volPath[] = new double[numberOfPaths];
		Arrays.fill(assetPath, Math.log(equityForwardStructure.getSpot()));
		Arrays.fill(volPath, sigma0);

		if(!saveMemory) {
			path[0][0] = assetPath;
			path[1][0] = volPath;
		}

		for(int i = 1; i < timeGrid.length; i++) {
			double deltaT = timeGrid[i] - timeGrid[i - 1];
			double sqrtDeltaT = Math.sqrt(deltaT);
			assert (sqrtDeltaT > 0) : "sqrt(delta) = 0!";
			double[][] increments = new double[numberOfPaths][2];
			// Fill Paths
			for(int j = 0; j < numberOfPaths; j++) {
				increments[j][0] = mersenneTwister.nextGaussian();

				// Vol
				double volPrev = volPath[j]; // Needed for asset
				double m = theta + (volPath[j] - theta) * Math.exp(-kappa * deltaT);
				double s2 = volPath[j] * epsilon * epsilon * Math.exp(-kappa * deltaT) / kappa * (1 - Math.exp(-kappa * deltaT))
						+ theta * epsilon * epsilon / (2 * kappa) * (1 - Math.exp(-kappa * deltaT)) * (1 - Math.exp(-kappa * deltaT));
				double psi = s2 / (m  * m);
				increments[j][1] = psi <= psiCritical ? mersenneTwister.nextGaussian() : mersenneTwister.nextDouble();
				if(psi <= psiCritical) {
					double b = Math.sqrt(2 / psi - 1 + Math.sqrt(2 / psi) * Math.sqrt(2 / psi - 1));
					double a = m / (1 + b * b);
					volPath[j] = a * (b + increments[j][1]) * (b + increments[j][1]);
				} else {
					double p = (psi - 1) / (psi + 1);
					double beta = (1 - p) / m;
					double bigPsiInverse = 0 <= increments[j][1] && increments[j][1] <= p ? 0 : 1 / beta * Math.log((1 - p) / (1 - increments[j][1]));
					volPath[j] = bigPsiInverse;
				}

				// Asset
				double k0 = - rho * kappa * theta / epsilon * deltaT;
				double k1 = gamma1 * deltaT * (kappa * rho / epsilon - 0.5) - rho / epsilon;
				double k2 = gamma2 * deltaT * (kappa * rho / epsilon - 0.5) + rho / epsilon;
				double k3 = gamma1 * deltaT + (1 - rho * rho);
				double k4 = gamma2 * deltaT * (1 - rho * rho);
				assetPath[j] = assetPath[j] + k0 + k1 * volPrev + k2 * volPath[j] + Math.sqrt(k3 * volPrev + k4 * volPath[j]) * increments[j][0];

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
				double avg = Arrays.stream(assetPathAtMaturities[currentMaturityIndex]).map(x -> Math.exp(x)).average().getAsDouble();
				for(int p = 0; p < numberOfPaths; p++) {
					assetPath[p] += Math.log(equityForwardStructure.getSpot() / avg);
					assetPathAtMaturities[currentMaturityIndex][p] = assetPath[p];
				}
				currentMaturityIndex++;
			}
		}
	}

}
