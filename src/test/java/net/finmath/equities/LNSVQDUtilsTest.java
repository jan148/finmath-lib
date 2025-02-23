package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LNSVQDUtilsTest {

	@Test
	public void matrixVectorMultTest() {
		Complex zero = Complex.ZERO;
		Complex one = Complex.ONE;

		Complex[][] matrix = new Complex[2][2];
		matrix[0][0] = one;
		matrix[0][1] = zero;
		matrix[1][0] = zero;
		matrix[1][1] = one;

		Complex[] vector = new Complex[2];
		vector[0] = new Complex(5, 0);
		vector[1] = new Complex(0, 5);

		Complex[] result = LNSVQDUtils.matrixVectorMult(matrix, vector);
		for(Complex num : result) {
			System.out.println(num.getReal() + "\t" + num.getImaginary());
		}
	}

	@Test
	public void addTimePointsToArrayTest() throws Exception {
		double[] array = {1., 2.45, 4.51, 8.7};
		int numPointsToInsert = 10;

		List<Double> arrayWithInsertedPoints = LNSVQDUtils.addTimePointsToArray(array, numPointsToInsert);

		System.out.println(arrayWithInsertedPoints.getClass().getSimpleName());
		System.out.println(arrayWithInsertedPoints);
	}

	@Test
	public void addMidPointsToListTest() throws Exception {
		ArrayList<Double> list = new ArrayList<>();
		list.add(1.);
		list.add(3.6);
		list.add(4.3);
		list.add(10.);

		List<Double> listWithMidPoints = LNSVQDUtils.addMidPointsToList(list);

		System.out.println(listWithMidPoints);
	}

	@Test
	public void printVDCNumbers() throws Exception {
		for(int j = 0; j < 100; j++) {
			System.out.println(LNSVQDUtils.modifiedVanDerCorput(j, 2));
		}
	}

	@Test
	public void sortListTest() throws Exception {
		int numberOfPoints = 22;
		int[] sortedInts = LNSVQDUtils.sortTimeIndices(numberOfPoints);
		for(int k = 0; k < numberOfPoints; k++) {
			System.out.println(k + " /  " + sortedInts[k]);
		}
	}

	// Creates a schedule for an index list 0, ..., numberOfPoints - 1
	@Test
	public void createSchedulingArrayTest() throws Exception {
		int numberOfPoints = 19;
		int[][] schedulingArray = LNSVQDUtils.createSchedulingArray(numberOfPoints);
		for(int k = 0; k < numberOfPoints; k++) {
			System.out.println(k + " \t  " + schedulingArray[k][0] + " \t  " + schedulingArray[k][1] + " \t  " + schedulingArray[k][2]);
		}
	}

	@Test
	public void checkBinString() {
		String x = LNSVQDUtils.getBinaryRepresentation(0.6763);
		System.out.println(x);
	}

	@Test
	public void confidenceIntervals() {
		StandardDeviation sd = new StandardDeviation(true);
		double[] arrayOld = {
				0.39216667, 0.39704215, 0.39423829, 0.40119638, 0.39436073,
				0.39663248, 0.39676409, 0.39810285, 0.39891461, 0.40348198
		};

		double[] array = {
				0.40107024, 0.39812148, 0.39653207, 0.39616844, 0.39676177,
				0.39665332, 0.39525346, 0.39836511, 0.39654689, 0.39959539
		}; // 400000 paths
		int numberOfPaths = 400000;
		double price = Arrays.stream(array).average().getAsDouble();
		double stdError = sd.evaluate(array);
		double lowerPrice = Math.max(price - 1.96 * stdError / Math.sqrt(array.length), 1E-10);
		double upperPrice = price + 1.96 * stdError / Math.sqrt(array.length);
					/*confidenceIntervalMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[2]);
					confidenceIntervalMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[3]);*/
		double[] confidenceInterval = new double[2];
		// confidenceInterval = LNSVQDUtils.getConfidenceInterval(array, 0.05);
		confidenceInterval[0] = lowerPrice;
		confidenceInterval[1] = upperPrice;
		System.out.println(price + "\t" + confidenceInterval[0] + "\t" + confidenceInterval[1]);
	}



}
