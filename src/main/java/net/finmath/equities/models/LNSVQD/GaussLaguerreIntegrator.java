package net.finmath.equities.models.LNSVQD;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialsUtils;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;

public class GaussLaguerreIntegrator {

	int n;

	PolynomialFunction laguerrePolynomial;
	Complex[] roots;

	public GaussLaguerreIntegrator(int n) {
		this.n = n;
		this.laguerrePolynomial = PolynomialsUtils.createLaguerrePolynomial(n);
	}

	/*double getRootsOfPolynomials() {
		LaguerreSolver laguerreSolver = new LaguerreSolver();
		roots = laguerreSolver.doSolve();
	}*/

	double integrater(UnivariateFunction univariateFunction) {
		double[] roots = new double[n];
		double integral = 0;
		for(int j = 0; j < n; j++) {
			double root = roots[j];
			double weight = root / Math.pow(n + 1, 2) ;
			double value = univariateFunction.value(root);
			integral += weight * value;
		}
		return integral;
		
	}

}
