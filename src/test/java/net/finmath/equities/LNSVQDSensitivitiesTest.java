package net.finmath.equities;

import net.finmath.equities.Simulation.Options.CliquetSimulationPricer;
import net.finmath.equities.models.LNSVQDUtils;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

public class LNSVQDSensitivitiesTest extends TestsSetupForLNSVQD {
	/**
	 * Stat utils
	 */
	StandardDeviation standardDeviation = new StandardDeviation();
	Random random = new Random();

	/*@Test
	public void testCliquetDeltaGamma() throws Exception {
		// Set the right case
		loadS25();

		int[] seeds = IntStream.range(0, 1).toArray();
		numberOfPaths = 100000;

		// Set Cliquet params
		double maturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double floorL = Double.NEGATIVE_INFINITY;
		double capL = 1.1;
		double floorG = 0;
		double capG = Double.POSITIVE_INFINITY;

		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();

		int minPercOfSpot = 50;
		int maxPercOfSpot = 150;
		double shiftSize = 1E-2; // 0.01

		double[] shiftPercs = IntStream.range(minPercOfSpot, maxPercOfSpot + 1).mapToDouble(n -> n / 100.).toArray();

		double[][] pricesLnsvqdMC = new double[seeds.length][shiftPercs.length];
		double[][] pricesLnsvqdQMC = new double[seeds.length][shiftPercs.length];
		double[][] pricesHestonMC = new double[seeds.length][shiftPercs.length];
		double[][] pricesHestonQMC = new double[seeds.length][shiftPercs.length];

		for(int s = 0; s < seeds.length; s++) {
			int seed = random.nextInt();

			// MC
			LNSVQDPathSimulatorMC pathSimulatorMC = new LNSVQDPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);;
			pathSimulatorMC.assetPathAtMaturities = new double[maturityGrid.length][numberOfPaths];
			Arrays.fill(pathSimulatorMC.assetPathAtMaturities[0], Math.log(spot0));

			// QMC
			LNSVQDPathSimulatorQMC pathSimulatorQMC = new LNSVQDPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, lnsvqdModelAnalyticalPricer, false);
			pathSimulatorQMC.assetPathAtMaturities = new double[maturityGrid.length][numberOfPaths];
			Arrays.fill(pathSimulatorQMC.assetPathAtMaturities[0], Math.log(spot0));

			double gamma1 = 1;
			double gamma2 = 0;

			// Heston MC
			HestonPathSimulatorMC pathSimulatorHestonMC = new HestonPathSimulatorMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
			pathSimulatorHestonMC.assetPathAtMaturities = new double[maturityGrid.length][numberOfPaths];
			Arrays.fill(pathSimulatorHestonMC.assetPathAtMaturities[0], Math.log(spot0));

			// Heston QMC
			HestonPathSimulatorQMC pathSimulatorHestonQMC = new HestonPathSimulatorQMC(valuationDate, disountCurve
					, equityForwardStructure, numberOfPaths, timeGrid, maturityGrid, selectedParamsHeston[0], selectedParamsHeston[1], selectedParamsHeston[2], selectedParamsHeston[3], selectedParamsHeston[4], gamma1, gamma2);
			pathSimulatorHestonQMC.assetPathAtMaturities = new double[maturityGrid.length][numberOfPaths];
			Arrays.fill(pathSimulatorHestonQMC.assetPathAtMaturities[0], Math.log(spot0));

			for(int k = 0; k < shiftPercs.length; k++) {
				double shift = shiftPercs[k];

				Boolean martingaleCorrection = Boolean.FALSE;
				int startingIndex = 93;
				double[] startingValueLNSVQD =  new double[]{spot0 * shift, selectedParamsLNSVQD[0]};
				double[] startingValueHeston =  new double[]{spot0 * shift, selectedParamsHeston[0]};

				// LNSVQD MC
				pathSimulatorMC.precalculatePaths(seed, true, startingIndex + 1, startingValueLNSVQD, martingaleCorrection);

				// LNSVQD QMC
				pathSimulatorQMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD, martingaleCorrection);

				// Heston MC
				pathSimulatorHestonMC.precalculatePaths(seed, true, startingIndex + 1, startingValueHeston, martingaleCorrection, );

				// Heston QMC
				pathSimulatorHestonQMC.precalculatePaths(seed, true, startingIndex, startingValueHeston, martingaleCorrection, );

				// Define pricers
				CliquetSimulationPricer simulPricerLNSVQDMC = new CliquetSimulationPricer<>(pathSimulatorMC);
				CliquetSimulationPricer simulPricerLNSVQDQMC = new CliquetSimulationPricer<>(pathSimulatorQMC);
				CliquetSimulationPricer simulPricerHestonMC = new CliquetSimulationPricer<>(pathSimulatorHestonMC);
				CliquetSimulationPricer simulPricerHestonQMC = new CliquetSimulationPricer<>(pathSimulatorHestonQMC);

				// MC LNSVQD
				double simulatedOptionPriceLNSVQDMC;
				try {
					simulatedOptionPriceLNSVQDMC = simulPricerLNSVQDMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
					simulatedOptionPriceLNSVQDMC = 1000000;
				}
				pricesLnsvqdMC[s][k] = simulatedOptionPriceLNSVQDMC;

				// QMC Heston
				double simulatedOptionPriceLNSVQDQMC;
				try {
					simulatedOptionPriceLNSVQDQMC = simulPricerLNSVQDQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
					simulatedOptionPriceLNSVQDQMC = 1000000;
				}
				pricesLnsvqdQMC[s][k] = simulatedOptionPriceLNSVQDQMC;

				// MC Heston
				double simulatedOptionPriceHestonMC;
				try {
					simulatedOptionPriceHestonMC = simulPricerHestonMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
					simulatedOptionPriceHestonMC = 1000000;
				}
				pricesHestonMC[s][k] = simulatedOptionPriceHestonMC;

				// QMC Heston
				double simulatedOptionPriceHestonQMC;
				try {
					simulatedOptionPriceHestonQMC = simulPricerHestonQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
					simulatedOptionPriceHestonQMC = 1000000;
				}
				pricesHestonQMC[s][k] = simulatedOptionPriceHestonQMC;
			}
			System.out.println("Done with seed: " + seed);
		}
		double[] pricesLnsvqdMCAvg = new double [shiftPercs.length];
		double[] pricesLnsvqdQMCAvg  = new double[shiftPercs.length];
		double[] pricesHestonMCAvg   = new double[shiftPercs.length];
		double[] pricesHestonQMCAvg  = new double[shiftPercs.length];
		for(int j = 0; j < shiftPercs.length; j++) {
			int index = j;
			pricesLnsvqdMCAvg[j] = IntStream.range(0, seeds.length).mapToDouble(i -> pricesLnsvqdMC[i][index]).average().getAsDouble();
			pricesLnsvqdQMCAvg[j] = IntStream.range(0, seeds.length).mapToDouble(i -> pricesLnsvqdQMC[i][index]).average().getAsDouble();
			pricesHestonMCAvg[j] = IntStream.range(0, seeds.length).mapToDouble(i -> pricesHestonMC[i][index]).average().getAsDouble();
			pricesHestonQMCAvg[j] = IntStream.range(0, seeds.length).mapToDouble(i -> pricesHestonQMC[i][index]).average().getAsDouble();
		}

		LNSVQDUtils.printArray(shiftPercs);

		double[][] firstAndSecondOrderLNMC = getFirstAndSecondOrderDerivatives(pricesLnsvqdMCAvg, shiftSize);
		double[][] firstAndSecondOrderLNQMC = getFirstAndSecondOrderDerivatives(pricesLnsvqdQMCAvg, shiftSize);
		double[][] firstAndSecondOrderHestonMC = getFirstAndSecondOrderDerivatives(pricesHestonMCAvg, shiftSize);
		double[][] firstAndSecondOrderHestonQMC = getFirstAndSecondOrderDerivatives(pricesHestonQMCAvg, shiftSize);

		for(int i = 0; i < 3; i++) {
			LNSVQDUtils.printArray(firstAndSecondOrderLNMC[i]);
			LNSVQDUtils.printArray(firstAndSecondOrderLNQMC[i]);
			LNSVQDUtils.printArray(firstAndSecondOrderHestonMC[i]);
			LNSVQDUtils.printArray(firstAndSecondOrderHestonQMC[i]);
		}
	}*/

	private double[][] getFirstAndSecondOrderDerivatives(double[] values, double shiftSize) {
		int numValues = values.length;

		double[][] output = new double[3][numValues];

		double[] firstDeriv = new double[numValues];
		for (int i = 1; i < numValues - 1; i++) {
			firstDeriv[i] = (values[i + 1] - values[i - 1]) / (2 * shiftSize);
		}

		double[] secondDeriv = new double[numValues];
		for (int i = 1; i < numValues - 1; i++) {
			secondDeriv[i] = (values[i + 1] - 2 * values[i] + values[i - 1]) / (shiftSize * shiftSize);
		}

		output[0] = values;
		output[1] = firstDeriv;
		output[2] = secondDeriv;

		return output;
	}

}