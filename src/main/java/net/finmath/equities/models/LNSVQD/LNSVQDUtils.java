package net.finmath.equities.models.LNSVQD;

import org.apache.commons.math3.complex.Complex;

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
		int rows = A.length; // Number of rows in the matrix
		int cols = A[0].length; // Number of columns in the matrix

		if (v.length != cols) {
			throw new IllegalArgumentException("Matrix column count must match vector size");
		}

		Complex[] result = new Complex[rows];
		for (int i = 0; i < rows; i++) {
			result[i] = new Complex(0, 0); // Start with zero
			for (int j = 0; j < cols; j++) {
				result[i] = result[i].add(A[i][j].multiply(v[j])); // Perform dot product
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


}
