package net.finmath.equities;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.SobolSequenceGenerator;
import org.junit.Test;

public class McQmcTest {

	@Test
	public void printMcQMC() {
		int numPts = 1000;
		SobolSequenceGenerator sobolSequenceGenerator = new SobolSequenceGenerator(2);
		MersenneTwister mersenneTwister = new MersenneTwister();
		sobolSequenceGenerator.nextVector();
		for(int i = 0; i < numPts; i++) {
			double[] sobolNum = sobolSequenceGenerator.nextVector();
			double[] mcNum = new double[2];
			mcNum[0] = mersenneTwister.nextDouble();
			mcNum[1] = mersenneTwister.nextDouble();
			System.out.println(mcNum[0] + "\t" + mcNum[1] + "\t" + sobolNum[0] + "\t" + sobolNum[1]);
		}
	}
}
