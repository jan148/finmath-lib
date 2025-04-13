package net.finmath.equities;

import net.finmath.equities.Simulation.HestonPathSimulator.HestonPathSimulatorMC;
import net.finmath.equities.Simulation.HestonPathSimulator.HestonPathSimulatorQMC;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorMC;
import net.finmath.equities.Simulation.LNSVQDPathSimulator.LNSVQDPathSimulatorQMC;
import net.finmath.equities.Simulation.Options.CliquetSimulationPricer;
import net.finmath.equities.models.LNSVQDUtils;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class LNSVQDSensitivitiesTest extends TestsSetupForLNSVQD {
	/**
	 * Stat utils
	 */
	StandardDeviation standardDeviation = new StandardDeviation();
	Random random = new Random();

	@Test
	public void testCliquetDelta() throws Exception {
		// Set the right case
		setDAXHestonMarchSetupSIM();

		int seed = 5555;
		numberOfPaths = 50000;

		// Set Cliquet params
		double maturity = strikeMatPairs.get(strikeMatPairs.size() - 1).getKey();
		double floorL = Double.NEGATIVE_INFINITY;
		double capL = 1.1;
		double floorG = 0;
		double capG = Double.POSITIVE_INFINITY;

		double[] timeGrid = LNSVQDUtils.addTimePointsToArray(maturityGrid,
						(int) (Math.round(maturity * 365.) * 1), 0, maturity, true)
				.stream().distinct().mapToDouble(Double::doubleValue).toArray();

		int minPercOfSpot = 30;
		int maxPercOfSpot = 300;
		double shiftSize = 1E-3;
		int incs = 10;

		double[] shiftBy = IntStream.range(0, maxPercOfSpot / incs).mapToDouble(n -> (minPercOfSpot + n * incs) / 100.).toArray();
		double[] shiftBy2 = Arrays.stream(shiftBy).map(x -> x + shiftSize).toArray();
		double[] shiftBy3 = Arrays.stream(shiftBy2).map(x -> x - 2 * shiftSize).toArray();
		double[] shiftPercs = DoubleStream.concat(
						DoubleStream.concat(Arrays.stream(shiftBy), Arrays.stream(shiftBy2)),
						Arrays.stream(shiftBy3)
				)
				.sorted()
				.toArray();

		double[] pricesLnsvqdMC = new double[shiftPercs.length];
		double[] pricesLnsvqdQMC = new double[shiftPercs.length];
		double[] pricesHestonMC = new double[shiftPercs.length];
		double[] pricesHestonQMC = new double[shiftPercs.length];

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

			int startingIndex = 94;
			double[] startingValueLNSVQD =  new double[]{spot0 * shift, selectedParamsLNSVQD[0]};
			double[] startingValueHeston =  new double[]{spot0 * shift, selectedParamsHeston[0]};

			// LNSVQD MC
			pathSimulatorMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD, Boolean.FALSE);

			// LNSVQD MC
			pathSimulatorQMC.precalculatePaths(seed, true, startingIndex, startingValueLNSVQD, Boolean.FALSE);

			// Heston MC
			pathSimulatorHestonMC.precalculatePaths(seed, true, startingIndex, startingValueHeston, Boolean.FALSE);

			// Heston QMC
			pathSimulatorHestonQMC.precalculatePaths(seed, true, startingIndex, startingValueHeston, Boolean.FALSE);

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
			pricesLnsvqdMC[k] = simulatedOptionPriceLNSVQDMC;

			// QMC Heston
			double simulatedOptionPriceLNSVQDQMC;
			try {
				simulatedOptionPriceLNSVQDQMC = simulPricerLNSVQDQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceLNSVQDQMC = 1000000;
			}
			pricesLnsvqdQMC[k] = simulatedOptionPriceLNSVQDQMC;

			// MC Heston
			double simulatedOptionPriceHestonMC;
			try {
				simulatedOptionPriceHestonMC = simulPricerHestonMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonMC = 1000000;
			}
			pricesHestonMC[k] = simulatedOptionPriceHestonMC;

			// QMC Heston
			double simulatedOptionPriceHestonQMC;
			try {
				simulatedOptionPriceHestonQMC = simulPricerHestonQMC.getCliquetPrice(maturity, floorL, capL, floorG, capG);
			} catch(AssertionError e) {
				System.err.println("Caught AssertionError: " + e.getMessage());
				simulatedOptionPriceHestonQMC = 1000000;
			}
			pricesHestonQMC[k] = simulatedOptionPriceHestonQMC;
		}

		LNSVQDUtils.printArray(shiftBy);

		double[][] firstAndSecondOrderLNMC = getFirstAndSecondOrderDerivatives(pricesLnsvqdMC, shiftSize);
		double[][] firstAndSecondOrderLNQMC = getFirstAndSecondOrderDerivatives(pricesLnsvqdQMC, shiftSize);
		double[][] firstAndSecondOrderHestonMC = getFirstAndSecondOrderDerivatives(pricesHestonMC, shiftSize);
		double[][] firstAndSecondOrderHestonQMC = getFirstAndSecondOrderDerivatives(pricesHestonQMC, shiftSize);

		for(int i = 0; i < 3; i++) {
			LNSVQDUtils.printArray(firstAndSecondOrderLNMC[i]);
			LNSVQDUtils.printArray(firstAndSecondOrderLNQMC[i]);
			LNSVQDUtils.printArray(firstAndSecondOrderHestonMC[i]);
			LNSVQDUtils.printArray(firstAndSecondOrderHestonQMC[i]);
		}
	}

	private double[][] getFirstAndSecondOrderDerivatives(double[] values, double shiftSize) {
		int numValues = values.length / 3;

		double[][] output = new double[3][numValues];
		double[] firstVals = new double[numValues];
		double[] secondVals = new double[numValues];
		double[] thirdVals  = new double[numValues];
		for (int i = 0; i < numValues; i++) {
			firstVals[i] = values[3 * i + 1];
			secondVals[i] = values[3 * i + 2];
			thirdVals[i] = values[3 * i];
		}

		double[] firstDeriv = new double[numValues];
		for (int i = 0; i < numValues; i++) {
			firstDeriv[i] = (secondVals[i] - firstVals[i]) / shiftSize;
		}

		double[] secondDeriv = new double[numValues];
		for (int i = 0; i < numValues; i++) {
			secondDeriv[i] = (secondVals[i] - 2 * firstVals[i] + thirdVals[i]) / (shiftSize * shiftSize);
		}

		output[0] = firstVals;
		output[1] = firstDeriv;
		output[2] = secondDeriv;

		return output;
	}

}