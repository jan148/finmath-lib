package net.finmath.equities;

import net.finmath.equities.models.LNSVQD.LNSVQDUtils;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

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
}
