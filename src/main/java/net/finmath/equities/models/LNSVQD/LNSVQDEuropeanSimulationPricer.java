package net.finmath.equities.models.LNSVQD;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class LNSVQDEuropeanSimulationPricer<T extends LNSVQDPathSimulator>{
	T pathSimulator;

	public LNSVQDEuropeanSimulationPricer(T pathSimulator) {
		this.pathSimulator = pathSimulator;
	}

	public double[] getPayoffsAtMaturity(double strike, double maturity, int callPutSign) throws Exception {
		List<Double> list = Arrays.stream(pathSimulator.maturities).boxed().collect(Collectors.toList());
		int matIndex = list.indexOf(maturity);
		if(matIndex == -1) {
			throw new Exception("Maturity not found!");
		}
		// TODO: Check next two lines
		double forwardFactor = pathSimulator.lnsvqdModel.equityForwardStructure.getForward(maturity) / pathSimulator.lnsvqdModel.spot0;
		double[] actualAssets = Arrays.stream(pathSimulator.assetPathAtMaturities[matIndex]).map(x -> Math.exp(x) * forwardFactor).toArray();
		double[] payoffsAtMaturity = Arrays.stream(actualAssets)
				.map(x -> Math.max(callPutSign * (x - strike), 0)).toArray();
		return payoffsAtMaturity;
	}

	public double getEuropeanPrice(double strike, double maturity, int callPutSign) throws Exception {
		double discountFactor = pathSimulator.lnsvqdModel.discountCurve.getDiscountFactor(maturity);
		double[] payoffsAtMaturity = getPayoffsAtMaturity(strike, maturity, callPutSign);
		// for(double payoff : payoffsAtMaturity) {assert(!Double.isNaN(payoff)) : "Nan encountered";}
		double expectationAtMaturity = Arrays.stream(payoffsAtMaturity).average().getAsDouble();
		double price = expectationAtMaturity * discountFactor;
		return price;
	}

	public double getEuropeanPriceAuto(double strike, double maturity) throws Exception {
		double forward = pathSimulator.lnsvqdModel.equityForwardStructure.getForward(maturity);
		int callPutSign;
		if(strike > forward) {
			callPutSign = 1;
		} else {
			callPutSign = -1;
		}
		return getEuropeanPrice(strike, maturity, callPutSign);
	}

}
