package net.finmath.equities.models.LNSVQD;

import java.util.Arrays;
import java.util.stream.IntStream;

public class CliquetSimulationPricer<T extends PathSimulator>{
	T pathSimulator;

	public CliquetSimulationPricer(T pathSimulator) {
		this.pathSimulator = pathSimulator;
	}

	public double[] getPayoffsAtMaturity(double maturity, double floorL, double capL, double floorG, double capG) throws Exception {
		double[] payoffsAtMaturity = new double[pathSimulator.numberOfPaths];

		for(int j = 0; j < pathSimulator.numberOfPaths; j++) {
			int pathIndex = j;
			double[] pathAtPayPoints = IntStream.range(0, pathSimulator.maturities.length)
					.mapToDouble(i ->
					{
						return Math.exp(pathSimulator.assetPathAtMaturities[i][pathIndex]);
					})
					.toArray();
			double payoff = 0;
			for(int m = 1; m < pathSimulator.maturities.length; m++) {
				payoff += Math.min(Math.max(floorL, pathAtPayPoints[m] / pathAtPayPoints[m - 1] - 1),  capL);
			}
			payoff = Math.min(Math.max(payoff, floorG), capG);
			payoffsAtMaturity[j] = payoff;
		}
		return payoffsAtMaturity;
	}

	public double getCliquetPrice(double maturity, double floorL, double capL, double floorG, double capG) throws Exception {
		double discountFactor = pathSimulator.discountCurve.getDiscountFactor(maturity);
		double[] payoffsAtMaturity = getPayoffsAtMaturity(maturity, floorL, capL, floorG, capG);
		// for(double payoff : payoffsAtMaturity) {assert(!Double.isNaN(payoff)) : "Nan encountered";}
		double expectationAtMaturity = Arrays.stream(payoffsAtMaturity).average().getAsDouble();
		double price = expectationAtMaturity * discountFactor;
		return price;
	}

}
