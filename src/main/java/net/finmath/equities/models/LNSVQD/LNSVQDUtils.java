package net.finmath.equities.models.LNSVQD;

import org.apache.commons.math3.complex.Complex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

public class LNSVQDUtils {
	/**
	 *
	 * Matrix operations
	 *
	 */
	public static Complex[][] matrixMult(Complex[][] A, Complex[][]  B) {
		int rowsA = A.length;
		int colsA = A[0].length;
		int rowsB = B.length;
		int colsB = B[0].length;

		if (colsA != rowsB) {
			throw new IllegalArgumentException("Matrix A's columns must equal Matrix B's rows");
		}

		Complex[][] result = new Complex[rowsA][colsB];
		for (int i = 0; i < rowsA; i++) {
			for (int j = 0; j < colsB; j++) {
				result[i][j] = new Complex(0, 0); // Initialize with 0 + 0i
				for (int k = 0; k < colsA; k++) {
					result[i][j] = result[i][j].add(A[i][k].multiply(B[k][j]));
				}
			}
		}
		return result;
	}

	public static Complex[] matrixVectorMult(Complex[][] A, Complex[]  v) {
		int rows = A.length;
		int cols = A[0].length;

		if (v.length != cols) {
			throw new IllegalArgumentException("Matrix column count must match vector size");
		}

		Complex[] result = new Complex[rows];
		for (int i = 0; i < rows; i++) {
			result[i] = new Complex(0, 0);
			for (int j = 0; j < cols; j++) {
				result[i] = result[i].add(A[i][j].multiply(v[j]));
			}
		}
		return result;
	}

	public static Complex scalarProduct(Complex[] v, Complex[]  w) {
		if (v.length != w.length) {
			throw new IllegalArgumentException("Vectors must be of the same length");
		}
		Complex result = new Complex(0, 0);
		for (int i = 0; i < v.length; i++) {
			result = result.add(v[i].multiply(w[i]));
			/*if(Double.isNaN(result.getReal()) || Double.isNaN(result.getImaginary())){
				System.out.println(v[i].getReal() + "\t" + v[i].getImaginary());
				System.out.println(w[i].getReal() + "\t" + w[i].getImaginary());
				throw new RuntimeException("STOP!");
			}*/
		}
		return result;
	}

	public static Complex[][] transposeMatrix(Complex[][] A) {
		int rows = A.length;
		int cols = A[0].length;

		Complex[][] transposed = new Complex[cols][rows];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				transposed[j][i] = A[i][j];
			}
		}
		return transposed;
	}

	/**
	 *
	 * Other
	 *
	 */
	public static double[] createTimeGrid(double t0, double tEnd, int steps) {
		double stepSize = (tEnd - t0) / steps;
		double[] timeGrid = new double[steps + 1];
		timeGrid[0] = t0;
		for(int j = 1; j <= steps; j++) {
			timeGrid[j] = timeGrid[j - 1] + stepSize;
		}
		return timeGrid;
	}
	public static void printArray(double[] array) {
		for(double element : array) {
			System.out.print(element + "\t");
		}
	}

	public static List<Double> addTimePointsToArray(double[] array, int numberOfInsertionPoints) throws Exception {
		double min = Arrays.stream(array).min().getAsDouble();
		double max = Arrays.stream(array).max().getAsDouble();
		double stepsize = (max - min) / (numberOfInsertionPoints + 1);

		double[] pointsToAdd = new double[numberOfInsertionPoints];
		pointsToAdd[0] = min + stepsize;
		for(int j = 1; j < numberOfInsertionPoints; j++) {
			pointsToAdd[j] = pointsToAdd[j - 1] + stepsize;
		}

		Double[] resultingArray = new Double[array.length + numberOfInsertionPoints];
		for(int k = 0; k < resultingArray.length; k++) {
			if(k < array.length) {
				resultingArray[k] = array[k];
			} else {
				resultingArray[k] = pointsToAdd[k - array.length];
			}
		}

		List<Double> list = Arrays.stream(resultingArray)
				.distinct()
				.sorted()
				.collect(Collectors.toList());

		if(list.size() != resultingArray.length) {
			throw new Exception("Adding points generated duplicates!");
		}

		return list;
	}

	public static List<Double> addMidPointsToList(List<Double> numbers) {
		List<Double> listWithMidpoints = new ArrayList<>();
		// Init
		listWithMidpoints.add(numbers.get(0));
		for (int i = 1; i < numbers.size(); i++) {
			double newNumber = numbers.get(i);
			double midPoint = numbers.get(i - 1) + ((newNumber - numbers.get(i - 1)) / 2.);
			listWithMidpoints.add(midPoint);
			listWithMidpoints.add(newNumber);
		}
		return listWithMidpoints;
	}


}
