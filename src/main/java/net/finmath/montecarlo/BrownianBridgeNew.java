package net.finmath.montecarlo;

import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import net.finmath.functions.NormalDistribution;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.time.TimeDiscretization;
import org.apache.commons.math3.random.SobolSequenceGenerator;

import java.io.IOException;
import java.util.*;

public class BrownianBridgeNew {
	private final ArrayList<Double> timeDiscretization;
	private transient double[][][] brownianIncrementsArray;
	private int[][] schedulingArray;
	private int numberOfPaths;

	public BrownianBridgeNew(final ArrayList<Double> timeDiscretization, int[][] schedulingArray, int numberOfPaths) {
		this.timeDiscretization = timeDiscretization;
		this.schedulingArray = schedulingArray;
		this.numberOfPaths = numberOfPaths;
	}

	public double[] getBrownianIncrementArr(final int timeIndex, final int factor, int seed) {
		// Thread safe lazy initialization
		if(brownianIncrementsArray == null) {
			doGenerateBrownianMotionNewFirstOverPaths(schedulingArray, seed);
		}

		/*
		 *  For performance reasons we return directly the stored data (no defensive copy).
		 *  We return an immutable object to ensure that the receiver does not alter the data.
		 */
		return brownianIncrementsArray[timeIndex][factor];
	}

	private void doGenerateBrownianMotionNewFirstOverPaths(int[][] schedulingArray, int seed) {
		Random random = new Random(seed); // Convert seed to long
		double scrambleNumber = random.nextDouble();
		assert (0. < scrambleNumber && scrambleNumber < 1.) : "ScrambleNumber is out of bounds!";
		// Initialize NGs; Need to check whether dim exceeds1 1000
		int rngDim = 2 * (schedulingArray.length - 1);
		SobolSequenceGenerator sobolSequenceGenerator = new SobolSequenceGenerator(Math.min(rngDim, 1000));
		sobolSequenceGenerator.nextVector();
		MersenneTwister mersenneTwister = new MersenneTwister(seed);

		if(brownianIncrementsArray != null) {
			return;
		}

		// Allocate memory
		brownianIncrementsArray = new double[timeDiscretization.size() - 1][getNumberOfFactors()][getNumberOfPaths()];

		double[][][] brownianMotionArr = new double[timeDiscretization.size()][getNumberOfFactors()][getNumberOfPaths()];
		double[][] initArray = new double[getNumberOfFactors()][getNumberOfPaths()];
		for(int i = 0; i < initArray.length; i++) {
			Arrays.fill(initArray[i], 0);  // Fill each row with 0s
		}
		brownianMotionArr[0] = initArray;

		for(int p = 0; p < getNumberOfPaths(); p++) {
			double[] vec = sobolSequenceGenerator.nextVector();
			double[] standardNormals = LNSVQDUtils.getStdNormalsFromUnifVec(vec, scrambleNumber);
			double[] standardAsset = new double[rngDim / 2];
			double[] standardVol = new double[rngDim / 2];
			for(int o = 0; o < standardNormals.length / 2; o++) {
				standardAsset[o] = standardNormals[2 * o];
				standardVol[o] = standardNormals[2 * o + 1];
			}

			if(rngDim > 1000) {
				for(int j = 500; j < rngDim / 2.; j++) {
					standardAsset[j] = NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());
					standardVol[j] = NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());
				}
			}

			for(int j = 1; j < brownianMotionArr.length; j++) {
				int timeIndex = schedulingArray[j][0];
				double time = timeDiscretization.get(timeIndex);

				int timeIndexEarlierTimeStep = schedulingArray[j][1];
				double timeEarlier = timeDiscretization.get(timeIndexEarlierTimeStep);

				int timeIndexLaterTimeStep = schedulingArray[j][2];
				double timeLater = timeDiscretization.get(timeIndexLaterTimeStep);

				double fac = (time - timeEarlier) / (timeLater - timeEarlier);

				for(int componentIndex = 0; componentIndex < 2; componentIndex++) {
					double[] standardNormal;
					if(componentIndex == 0) {
						standardNormal = standardAsset;
					} else {
						standardNormal = standardVol;
					}

					if(j == 1) {
						brownianMotionArr[timeIndex][componentIndex][p] =
								brownianMotionArr[timeIndexEarlierTimeStep][componentIndex][p] + standardNormal[j - 1] * Math.sqrt(time);
					} else {
						brownianMotionArr[timeIndex][componentIndex][p] =
								brownianMotionArr[timeIndexEarlierTimeStep][componentIndex][p]
										+ ((brownianMotionArr[timeIndexLaterTimeStep][componentIndex][p] - brownianMotionArr[timeIndexEarlierTimeStep][componentIndex][p]) * fac)
										+ (standardNormal[j - 1] * Math.sqrt(fac * (timeLater - time)));
					}
				}
			}
		}

		double[] brownianMotionLastTime; //  = brownianMotion[j + 1][i].getRealizations();
		double[] brownianMotionCurrentTime;
		for(int j = 0; j < timeDiscretization.size() - 1; j++) {
			for(int i = 0; i < getNumberOfFactors(); i++) {
				brownianMotionLastTime = brownianMotionArr[j][i];
				brownianMotionCurrentTime = brownianMotionArr[j + 1][i];
				for(int p = 0; p < getNumberOfPaths(); p++) {
					brownianIncrementsArray[j][i][p] = brownianMotionCurrentTime[p] - brownianMotionLastTime[p];
				}
			}
		}
	}

	public int getNumberOfFactors() {
		return 2;
	}

	public int getNumberOfPaths() {
		return numberOfPaths;
	}

}
