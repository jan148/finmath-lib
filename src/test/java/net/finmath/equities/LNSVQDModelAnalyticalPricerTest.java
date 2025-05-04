package net.finmath.equities;

import net.finmath.equities.models.LNSVQDUtils;
import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.stream.IntStream;

import static org.junit.Assert.assertEquals;

/**
 * This test
 * 1. compares the semi-analytical LNSVQD call option price to the BS call option price and
 * 2. compares the semi-analytical LNSVQD call option price to the simulated LNSVQD call option price
 */
class LNSVQDModelAnalyticalPricerTest extends TestsSetupForLNSVQD{
	/**
	 * Tolerance level
	 */
	private final double delta = 10E-3;


	/**
	 * ***************************************************+
	 * 1. Comparison semi-analytical LNSVQD call <-> BS call
	 * ***************************************************+
	 */
	@Test
	public void printAj() {
		int index = 3;
		double ttm = 1;
		double y = 2;
		int numStepsForODEIntegration = (int) (ttm * 365 * lnsvqdModelAnalyticalPricer.getNumStepsForODEIntegrationPerYear());
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, numStepsForODEIntegration);
		final Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
		Complex[][] solutionPath = lnsvqdModelAnalyticalPricer.getSolutionPathForODESystem(timeGrid, charFuncArgs);
		System.out.println("T \t Real part \t Imaginary part");
		for(int i = 0; i < solutionPath.length; i++) {
			double t = timeGrid[i];
			Complex value = solutionPath[i][index];
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
		}
	}

	@Test
	public void printE2() {
		loadBitcoin();

		double ttm = 1;
		double y = -1.085746;
		int numStepsForODEIntegration = (int) (ttm * lnsvqdModelAnalyticalPricer.getNumStepsForODEIntegrationPerYear());
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, numStepsForODEIntegration);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximation(t, charFuncArgs)
					.multiply(charFuncArgs[0].multiply(lnsvqdModelAnalyticalPricer.getX0())
							.add(charFuncArgs[1].multiply(lnsvqdModelAnalyticalPricer.getI0())).exp());
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
		}
	}

	@Test
	public void printMGF() {
		double ttm = 0;
		double y = 2;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, 100);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximation(t, charFuncArgs);
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
		}
	}

	// Next test work only if volatility parameters are 0, i.e. model is a GBM
	@Test
	void calculateExponentialAffineApproximation() {
		double maturity = 1;
		if(selectedParamsLNSVQD[1] != 0 || selectedParamsLNSVQD[2] != 0 || selectedParamsLNSVQD[3] != 0 || selectedParamsLNSVQD[4] != 0 || selectedParamsLNSVQD[5] != 0) {
			return;
		}
		double y = 1.;
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};
		Complex exponentialAffineApproximationAnalyticalValue = (charFuncArgs[0].multiply(-lnsvqdModelAnalyticalPricer.getX0())
				.add(charFuncArgs[0].pow(2).add(charFuncArgs[0]).subtract(charFuncArgs[1].multiply(2)).multiply(0.5 * maturity)
						.multiply(Math.pow(lnsvqdModelAnalyticalPricer.getY0(), 2)))).exp();
		Complex exponentialAffineApproximationOdeValue = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximation(maturity, charFuncArgs);
		System.out.println("Analytical exponential-affine approximation value at " + maturity + ": " + exponentialAffineApproximationAnalyticalValue);
		System.out.println("ODE-based exponential-affine approximation value at " + maturity + ": " + exponentialAffineApproximationOdeValue);

		assertEquals(exponentialAffineApproximationAnalyticalValue.getReal(), exponentialAffineApproximationOdeValue.getReal(), delta);
		assertEquals(exponentialAffineApproximationAnalyticalValue.getImaginary(), exponentialAffineApproximationOdeValue.getImaginary(), delta);
	}

	/**
	 * ***************************************************+
	 * SECTION 3: TEST CALL PRICE CALCULATION WITH PRECALCULATION OF E2
	 * ***************************************************+
	 */
	@Test
	public void calculateExponentialAffineApproximationFullPathTest() {
		double ttm = 1;
		double y = 2.;
		int numStepsForODEIntegration = (int) (ttm * 365 * lnsvqdModelAnalyticalPricer.getNumStepsForODEIntegrationPerYear());
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0, ttm, numStepsForODEIntegration);
		Complex[] charFuncArgs = new Complex[]{new Complex(-0.5, y), Complex.ZERO, Complex.ZERO};

		System.out.println("TTM \t Real part \t Imaginary part");
		Complex[] values = lnsvqdModelAnalyticalPricer.calculateExponentialAffineApproximationFullPath(timeGrid, charFuncArgs);
		for(int i = 0; i < timeGrid.length; i++) {
			double t = timeGrid[i];
			Complex value = values[i].multiply(charFuncArgs[0].multiply(lnsvqdModelAnalyticalPricer.getX0()).add(charFuncArgs[1].multiply(lnsvqdModelAnalyticalPricer.getI0())).exp());
			double realPart = value.getReal();
			double imaginaryPart = value.getImaginary();
			System.out.println(t + "\t" + realPart + "\t" + imaginaryPart);
		}
	}

	/**
	 * ***************************************************+
	 * SECTION 5: TEST DEPENDENCE ON HYPERPARAMS
	 * ***************************************************+
	 */
	@Test
	public void printPricesWrtPricerParams() throws Exception {
		loadS24();
		int[] vals = IntStream.range(1, 11).map(i -> 40 * i).toArray();
		double[][] impliedVolsPerStepCounter = new double[vals.length][strikeMatPairs.size()];
		for(int j = 0; j < vals.length; j++) {
			lnsvqdModelAnalyticalPricer.setUpperBoundForIntegration(vals[j]); // lnsvqdModelAnalyticalPricer.numStepsForODEIntegrationPerYear = vals[j];
			impliedVolsPerStepCounter[j] = lnsvqdModelAnalyticalPricer.getImpliedVolSurfaceFromStrikeMatList(strikeMatPairs);
			LNSVQDUtils.printArray(impliedVolsPerStepCounter[j]);
		}
		for(int i = 0; i < strikeMatPairs.size(); i++) {
			double[] res = new double[vals.length];
			for(int k = 0; k < vals.length; k++) {
				res[k] = impliedVolsPerStepCounter[k][i];
			}
			double min =  Arrays.stream(res).min().getAsDouble();
			double max =  Arrays.stream(res).max().getAsDouble();
			double range = max - min;

			assertEquals(0., range, 1E-4);
		}
	}

}
