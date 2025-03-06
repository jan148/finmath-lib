package net.finmath.integration;

import net.finmath.equities.models.LNSVQDUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegratorFactory;
import org.junit.Test;

import java.util.ArrayList;

public class IterativeLegendreGaussIntergratorTest {
	@Test
	public void printLegendrePointsAndIntegrationValue() {
		// Create an integrator instance
		// Parameters: number of points, relative accuracy, absolute accuracy, min and max iterations
		int numberOfPoints = 15; // Number of integration points for Legendre-Gauss quadrature
		double relativeAccuracy = 1.0e-6;
		double absoluteAccuracy = 1.0e-9;
		int minIterations = 3;
		int maxIterations = 100;

		UnivariateFunction function = x -> (1 / (x * x + 0.25));

		IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(
				numberOfPoints, relativeAccuracy, absoluteAccuracy, minIterations, maxIterations
		);

		// Define the integration interval [a, b]
		double lowerBound = 0;
		double upperBound = 20;

		// Print Legendre points
		// We need the factory to get the Legender points
		final GaussIntegratorFactory FACTORY = new GaussIntegratorFactory();

		// The GL-implementation computes the solution on equally-sized subintervals
		ArrayList<Double> legendrePoints = new ArrayList<>();
		double[] subintervalPoints = LNSVQDUtils.createTimeGrid(lowerBound, upperBound, numberOfPoints);
		double stepSize = subintervalPoints[1] - subintervalPoints[0];
		for(int i = 0; i < numberOfPoints - 1; i++) {
			double a = subintervalPoints[0];
			double b = subintervalPoints[1];
			double c = b - a;
			final GaussIntegrator g = FACTORY.legendreHighPrecision(numberOfPoints, a, b);
			legendrePoints.add(g.getPoint(i));
		}

		System.out.println(legendrePoints.toString());

		// Perform the integration
		double result = 0;
		try {
			result = integrator.integrate(100000000, function, lowerBound, upperBound);
			System.out.println("Integral result: " + result);
		} catch (Exception e) {
			System.err.println("Integration failed: " + e.getMessage());
		}
	}


}
