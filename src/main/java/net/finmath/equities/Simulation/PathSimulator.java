package net.finmath.equities.Simulation;

import net.finmath.equities.marketdata.YieldCurve;
import net.finmath.equities.models.EquityForwardStructure;

import java.time.LocalDate;

public abstract class PathSimulator {
	public LocalDate spotDate;
	public YieldCurve discountCurve;
	public EquityForwardStructure equityForwardStructure;
	public int numberOfPaths;
	public double[] timeGrid;
	public double[] maturities;
	public double path[][][];
	public double assetPathAtMaturities[][];
	public abstract void precalculatePaths(int seed, Boolean saveMemory);

}
