package net.finmath.equities;

import net.finmath.equities.Simulation.PathSimulator;
import net.finmath.equities.models.LNSVQDUtils;
import org.junit.Test;

public class LNSVQDPaths extends TestsSetupForLNSVQD {

	@Test
	public void rollOutPathsLNSVQD() throws Exception {
		loadS24();
		numberOfPaths = 10;
		int seed = 7205;
		double maxMaturity = 10; // strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] startingValueLNSVQD = new double[]{spot0, selectedParamsLNSVQD[0]};
		double[] extended = new double[maturityGrid.length + 1];
		System.arraycopy(maturityGrid, 0, extended, 0, maturityGrid.length);
		extended[extended.length - 1] = maxMaturity;
		maturityGrid = extended;
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		PathSimulator pathSimulatorMC = new PathSimulator(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid);
		pathSimulatorMC.precalculatePaths(seed, false, 1, startingValueLNSVQD, Boolean.TRUE, "LNSVQD", "MC", lnsvqdModelDescriptor, null);
		for(int t = 0; t < timeGrid.length; t++) {
			System.out.print(timeGrid[t] + "\t");
			for(int k = 0; k < numberOfPaths; k++) {
				System.out.print(pathSimulatorMC.path[1][t][k] + "\t");
			}
			System.out.println();
		}
	}

	@Test
	public void rollOutPathsHeston() throws Exception {
		loadS24();
		numberOfPaths = 10;
		int seed = 7205;
		double maxMaturity = 10; // strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] startingValueHeston = new double[]{spot0, selectedParamsHeston[0]};
		double[] extended = new double[maturityGrid.length + 1];
		System.arraycopy(maturityGrid, 0, extended, 0, maturityGrid.length);
		extended[extended.length - 1] = maxMaturity;
		maturityGrid = extended;
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		// Heston MC
		double gamma1 = 1;
		double gamma2 = 0;
		PathSimulator pathSimulatorMC = new PathSimulator(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid);
		pathSimulatorMC.precalculatePaths(seed, false, 1, startingValueHeston, Boolean.TRUE, "Heston", "MC", null, hestonModelDescriptor);
		for(int t = 0; t < timeGrid.length; t++) {
			System.out.print(timeGrid[t] + "\t");
			for(int k = 0; k < numberOfPaths; k++) {
				System.out.print(Math.sqrt(pathSimulatorMC.path[1][t][k]) + "\t");
			}
			System.out.println();
		}
	}
}
