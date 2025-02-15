package net.finmath.equities.models.LNSVQD;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Pair;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
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
		// TODO: A cheat, change later; problem is precision
		timeGrid[steps] = tEnd;
		return timeGrid;
	}

	public static void printArray(double[] array) {
		System.out.print("\n");
		for(double element : array) {
			System.out.print(element + "\t");
		}
		System.out.print("\n");
	}

	public static void printPricesFromMaturityStrikeGrid(double[] maturities, double[] strikes, double[] prices) {
		for(int s = 0; s < strikes.length; s++) {
			System.out.print("\t" + strikes[s]);
		}
		System.out.println();
		for(int m = 0; m < maturities.length; m++) {
			double maturity = maturities[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(prices[m * strikes.length + s] + "\t");
			}
			System.out.print("\n");
		}
	}

	// TODO
	public static List<Double> addTimePointsToArray(double[] array, int numberOfInsertionPoints) throws Exception {
		double from = Arrays.stream(array).min().getAsDouble();
		double to = Arrays.stream(array).max().getAsDouble();
		return addTimePointsToArray(array, numberOfInsertionPoints, from, to, false);
	}

	// fromToIncl = from / to inclusive?
	public static List<Double> addTimePointsToArray(double[] array, int numberOfInsertionPoints, double from, double to, Boolean fromToIncl) throws Exception {
		double min = from;
		double max = to;
		double stepsize = fromToIncl ? (max - min) / numberOfInsertionPoints : (max - min) / (numberOfInsertionPoints + 1);

		double[] pointsToAdd = new double[numberOfInsertionPoints];
		pointsToAdd[0] = fromToIncl ? min : min + stepsize;
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

		if(array.length + numberOfInsertionPoints != list.size()) {
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

	public static Complex[] convertDoubleArrayToComplexArray(double[] array) {
		Complex[] complexes = new Complex[array.length];
		for (int j = 0; j < array.length; j++) {
			complexes[j] = new Complex(array[j], 0.);
		}
		return complexes;
	}

	public static double[] convertComplexArrayToDoubleWithReal(Complex[] array) {
		double[] reals = new double[array.length];
		for (int j = 0; j < array.length; j++) {
			reals[j] = array[j].getReal();
		}
		return reals;
	}

	public static List<Pair<Double, Double>> create2dMesh(double[] arr1, double[] arr2) {
		List<Pair<Double, Double>> mesh= new ArrayList<>();
		for (int i = 0; i < arr1.length; i++) {
			for (int j = 0; j < arr2.length; j++) {
				mesh.add(new Pair<Double, Double>(arr1[i], arr2[j]));
			}
		}
		return mesh;
	}

	public static LocalDate[] createLocalDateList(String[] dates) {
		LocalDate[] datesArray = new LocalDate[dates.length];
		for (int i = 0; i < dates.length; i++) {
			datesArray[i] = LocalDate.parse(dates[i]);
		}
		return datesArray;
	}
}
