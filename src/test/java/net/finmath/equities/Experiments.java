package net.finmath.equities;

import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

public class Experiments {

	@Test
	public void calcComplexValues() {
		Complex x = new Complex(1E-30, 0);
		System.out.println(x.multiply(1E-300));
	}
}
