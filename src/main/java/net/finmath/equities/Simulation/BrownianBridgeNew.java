package net.finmath.equities.Simulation;

import net.finmath.functions.NormalDistribution;
import net.finmath.randomnumbers.MersenneTwister;

import java.util.*;

/**
 * Construction of a Brownian bridg; The order of construction of the Brownian observations are
 * determined by the scheduling array.
 *
 * @author Jan Berger
 */
public class BrownianBridgeNew {
	private final ArrayList<Double> timeDiscretization;
	private int[][] schedulingArray;
	private int numberOfPaths;

	public BrownianBridgeNew(final ArrayList<Double> timeDiscretization, int[][] schedulingArray, int numberOfPaths) {
		this.timeDiscretization = timeDiscretization;
		if(schedulingArray.length < 2) {
			throw new Error("schedulingArray must contain at least two points!");
		}
		this.schedulingArray = schedulingArray;
		this.numberOfPaths = numberOfPaths;
	}

	public double[][] generateBrownianIncrementsOnePath(double[] normals, MersenneTwister mersenneTwister) {
		double[][] brownianMotionArr = new double[timeDiscretization.size()][2];
		for(int i = 0; i < schedulingArray.length; i++) {
			brownianMotionArr[i][0] = 0;
			brownianMotionArr[i][1] = 0;
		}
		double[][] brownianIncs = new double[timeDiscretization.size() - 1][2];

		int rngDim = 2 * (schedulingArray.length - 1);

		double[] standardAsset = new double[rngDim / 2];
		double[] standardVol = new double[rngDim / 2];
		for(int o = 0; o < normals.length / 2; o++) {
			standardAsset[o] = normals[2 * o];
			standardVol[o] = normals[2 * o + 1];
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
					brownianMotionArr[timeIndex][componentIndex] =
							brownianMotionArr[timeIndexEarlierTimeStep][componentIndex] + standardNormal[j - 1] * Math.sqrt(time - timeEarlier);
				} else {
					brownianMotionArr[timeIndex][componentIndex] =
							brownianMotionArr[timeIndexEarlierTimeStep][componentIndex]
									+ ((brownianMotionArr[timeIndexLaterTimeStep][componentIndex] - brownianMotionArr[timeIndexEarlierTimeStep][componentIndex]) * fac)
									+ (standardNormal[j - 1] * Math.sqrt(fac * (timeLater - time)));
				}
			}
		}

		double brownianMotionLastTime; //  = brownianMotion[j + 1][i].getRealizations();
		double brownianMotionCurrentTime;
		for(int j = 0; j < timeDiscretization.size() - 1; j++) {
			for(int i = 0; i < getNumberOfFactors(); i++) {
				brownianMotionLastTime = brownianMotionArr[j][i];
				brownianMotionCurrentTime = brownianMotionArr[j + 1][i];
				brownianIncs[j][i] = brownianMotionCurrentTime - brownianMotionLastTime;
			}
		}

		return brownianIncs;
	}

	public int getNumberOfFactors() {
		return 2;
	}

	public int getNumberOfPaths() {
		return numberOfPaths;
	}

}
