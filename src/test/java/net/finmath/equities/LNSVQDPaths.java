package net.finmath.equities;

import net.finmath.equities.Simulation.HestonPathSimulator.HestonPathSimulatorMC;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorMC;
import net.finmath.equities.models.LNSVQDUtils;
import org.junit.Test;

public class LNSVQDPaths extends TestsSetupForLNSVQD {

	@Test
	public void rollOutPathsLNSVQD() throws Exception {
		setDAXHestonMarchSetupSIM();
		int seed = 7205;
		double maxMaturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
		pathSimulatorMC.precalculatePaths(seed, false);
		for(int t = 0; t < timeGrid.length; t++) {
			for(int k = 0; k < numberOfPaths; k++) {
				System.out.print(pathSimulatorMC.path[1][t][k] + "\t");
			}
			System.out.println();
		}
	}

	@Test
	public void rollOutPathsHeston() throws Exception {
		setDAXHestonMarchSetupSIM();
		int seed = 7205;
		double maxMaturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maxMaturity * 365.) * 1), 0, maxMaturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();
		// Heston MC
		double gamma1 = 1;
		double gamma2 = 0;
		HestonPathSimulatorMC pathSimulatorHestonMC = new HestonPathSimulatorMC(valuationDate, disountCurve
				, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
		// sw.reset();
		// sw.start();
		pathSimulatorHestonMC.precalculatePaths(seed, false);
		for(int t = 0; t < timeGrid.length; t++) {
			for(int k = 0; k < numberOfPaths; k++) {
				System.out.print(pathSimulatorHestonMC.path[1][t][k] + "\t");
			}
			System.out.println();
		}
	}
}
