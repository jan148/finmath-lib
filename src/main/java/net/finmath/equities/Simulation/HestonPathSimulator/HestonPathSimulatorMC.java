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

			// Pre-chaching
			// Asset
			double A = rho / epsilon * (1 + kappa * gamma2 * deltaT) - 0.5 * gamma2 * deltaT * rho * rho;
			double k0 = - rho * kappa * theta / epsilon * deltaT;
			double k0Star;
			double k1 = gamma1 * deltaT * (kappa * rho / epsilon - 0.5) - rho / epsilon;
			double k2 = gamma2 * deltaT * (kappa * rho / epsilon - 0.5) + rho / epsilon;
			double k3 = gamma1 * deltaT + (1 - rho * rho);
			double k4 = gamma2 * deltaT * (1 - rho * rho);

			// Fill Paths
			for(int j = 0; j < numberOfPaths; j++) {
				increments[j][0] = mersenneTwister.nextGaussian();

				// Vol
				double volPrev = volPath[j]; // Needed for asset
				double m = theta + (volPrev - theta) * Math.exp(-kappa * deltaT);
				double s2 = volPrev * epsilon * epsilon * Math.exp(-kappa * deltaT) / kappa * (1 - Math.exp(-kappa * deltaT))
						+ theta * epsilon * epsilon / (2 * kappa) * (1 - Math.exp(-kappa * deltaT)) * (1 - Math.exp(-kappa * deltaT));
				double psi = s2 / (m  * m);
				if(psi <= psiCritical) {
					increments[j][1] = mersenneTwister.nextGaussian();
					double b = Math.sqrt(2 / psi - 1 + Math.sqrt(2 / psi) * Math.sqrt(2 / psi - 1));
					double a = m / (1 + b * b);
					volPath[j] = a * (b + increments[j][1]) * (b + increments[j][1]);
					k0Star = -A * b * b * a / (1 - 2 * A * a) + 0.5 * Math.log(1 - 2 * A * a) - (k1 + 0.5 * k3) * volPrev;
				} else {
					increments[j][1] = mersenneTwister.nextDouble();
					double p = (psi - 1) / (psi + 1);
					double beta = (1 - p) / m;
					// System.out.println((1 - p) / beta + "\t" + m);
					// System.out.println((1 - p * p) / (beta * beta) + "\t" + s2);
					double bigPsiInverse = 0 <= increments[j][1] && increments[j][1] <= p ? 0 : 1 / beta * Math.log((1 - p) / (1 - increments[j][1]));
					volPath[j] = bigPsiInverse;
					k0Star = -Math.log(p + beta * (1 - p) / (beta - A)) - (k1 + 0.5 * k3) * volPrev;
				}

				// Asset
				assetPath[j] = assetPath[j] + k0Star + k1 * volPrev + k2 * volPath[j] + Math.sqrt(k3 * volPrev + k4 * volPath[j]) * increments[j][0];

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
			// Increment maturity index
			if(maturities[currentMaturityIndex] == timeGrid[i]) {
				currentMaturityIndex++;
			}
		}
	}

}
