package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import java.util.ArrayList;
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
		int numberOfPoints = 22;
		int[][] schedulingArray = LNSVQDUtils.createSchedulingArray(numberOfPoints);
		for(int k = 0; k < numberOfPoints; k++) {
			System.out.println(k + " /  " + schedulingArray[k][0] + " /  " + schedulingArray[k][1] + " /  " + schedulingArray[k][2]);
		}
	}


}
