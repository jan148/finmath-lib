package net.finmath.equities.models.LNSVQD;

import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;

import java.time.LocalDate;

public abstract class PathSimulator {
	LocalDate spotDate;
	YieldCurve discountCurve;
	EquityForwardStructure equityForwardStructure;
	int numberOfPaths;
	double[] timeGrid;
	double[] maturities;
	double path[][][];
	double assetPathAtMaturities[][];
	public abstract void precalculatePaths(int seed, Boolean saveMemory);

}
