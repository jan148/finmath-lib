package net.finmath.equities.models.LNSVQD;

import org.apache.commons.math3.analysis.UnivariateFunction;

public abstract class LNSVQDPathSimulator {
	LNSVQDModel lnsvqdModel;
	int numberOfPaths;
	double[] timeGrid;
	double[] maturities;
	double path[][][];
	double assetPathAtMaturities[][];
	boolean isBackwardEuler = false;
	UnivariateFunction zeta = x -> Math.exp(-x) * lnsvqdModel.getKappa1() * lnsvqdModel.getTheta() - Math.exp(x) * lnsvqdModel.getKappa2()
			- lnsvqdModel.getKappa1() + lnsvqdModel.getKappa2() * lnsvqdModel.getTheta() - 0.5 * lnsvqdModel.getTotalInstVar();

	public LNSVQDPathSimulator(LNSVQDModel lnsvqdModel, int numberOfPaths, double[] timeGrid, double[] maturities, Boolean isBackwardEuler) {
		this.lnsvqdModel = lnsvqdModel;
		this.numberOfPaths = numberOfPaths;
		this.timeGrid = timeGrid;
		this.maturities = maturities;
		// Component, time, paths
		this.isBackwardEuler = isBackwardEuler;
	}

	public abstract void precalculatePaths(int seed, Boolean saveMemory);

}
