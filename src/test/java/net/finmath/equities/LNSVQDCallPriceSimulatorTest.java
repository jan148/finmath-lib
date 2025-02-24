package net.finmath.equities;

import net.finmath.equities.marketdata.VolatilityPoint;
import net.finmath.equities.models.Black76Model;
import net.finmath.equities.models.LNSVQD.*;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.*;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloLNSVQDModel;
import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.montecarlo.process.LNSVQDDiscretizationScheme;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import java.time.Duration;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class LNSVQDCallPriceSimulatorTest extends TestsSetupForLNSVQD {
	/**
	 * Stat utils
	 */
	StandardDeviation standardDeviation = new StandardDeviation();
	Random random = new Random(3172);

	@Test
	public void getCallPrices() throws Exception {
		// Get option values
		double[] strikes = LNSVQDUtils.createTimeGrid(0.6 * spot0, 1.4 * spot0, 4);
		List<Pair<Double, Double>> strikeMaturityPairs = LNSVQDUtils.create2dMesh(maturityGrid, strikes);
		//double[] relativeErrors = new double[maturityGrid.length];

		double[][] pricesBS = new double[maturityGrid.length][strikes.length];
		double[][] pricesMC = new double[maturityGrid.length][strikes.length];
		double[][] pricesQMC = new double[maturityGrid.length][strikes.length];
		double[][] pricesMCFinmath = new double[maturityGrid.length][strikes.length];
		// Get analytical prices
		double[] pricesAnalytical = lnsvqdModelAnalyticalPricer.getEuropeanOptionPrices(strikeMaturityPairs, true);

		double[][] stdErrorsMc = new double[maturityGrid.length][strikes.length];
		double[][] stdErrorsMcFinmath = new double[maturityGrid.length][strikes.length];
		double[][] stdErrorsQMc = new double[maturityGrid.length][strikes.length];

		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			double discountFactor = equityForwardStructure.getRepoCurve().getDiscountFactor(maturity);
			double forward = equityForwardStructure.getForward(maturity);
			double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
					maturity, (int) Math.round(maturity * 365.));
			for(int s = 0; s < strikes.length; s++) {
				double strike = strikes[s];

				// BS value
				pricesBS[m][s] = Black76Model.optionPrice(forward, strike, maturity, selectedParams[0], true, discountFactor);

				List<Integer> seeds = random.ints(10).boxed().collect(Collectors.toList());
				double[] prices = new double[seeds.size()];
				double[] pricesFm = new double[seeds.size()];
				double[] pricesQ = new double[seeds.size()];

				for(int seed : seeds) {
					// Normal MC
					LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdCallPriceSimulator.precalculatePaths(seed);
					double simulatedOptionPrice = lnsvqdCallPriceSimulator.getCallPrice(strike, maturity, 1);
					prices[seeds.indexOf(seed)] = simulatedOptionPrice;

					// QMC
					LNSVQDPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
					double simulatedOptionPriceQMC = lnsvqdPriceSimulatorQMC.getCallPrice(strike, maturity, 1);
					pricesQ[seeds.indexOf(seed)] = simulatedOptionPriceQMC;

					lnsvqdPriceSimulatorQMC = null;

					// Price for finmath implementation
					// BrownianMotionFromMersenneRandomNumbers brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
					BrownianMotionFromMersenneRandomNumbers brownianBridge = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, 2, numberOfPaths, seed, randomVariableFactory);
					LNSVQDDiscretizationScheme lnsvqdDiscretizationScheme = new LNSVQDDiscretizationScheme(lnsvqdModelAnalyticalPricer, brownianBridge);
					MonteCarloLNSVQDModel monteCarloLNSVQDModel = new MonteCarloLNSVQDModel(lnsvqdDiscretizationScheme, seed);
					EuropeanOption europeanOption = new EuropeanOption(maturity, strike, 1, 0);
					pricesFm[seeds.indexOf(seed)] = 123456789; // europeanOption.getValue(monteCarloLNSVQDModel);
				}

				double averagePrice = Arrays.stream(prices).average().getAsDouble();
				double varMC = Arrays.stream(prices).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
				double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());

				double averagePriceFm = Arrays.stream(pricesFm).average().getAsDouble();
				double varMcFm = Arrays.stream(pricesFm).map(x -> Math.pow(x - averagePriceFm, 2)).sum() / (seeds.size() - 1);
				double stdErrMcFm = Math.sqrt(varMcFm) / Math.sqrt(seeds.size());

				double averagePriceQMC = Arrays.stream(pricesQ).average().getAsDouble();
				double varQMC = Arrays.stream(pricesQ).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
				double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());

				pricesMC[m][s] = averagePrice;
				pricesMCFinmath[m][s] = averagePriceFm;
				pricesQMC[m][s] = averagePriceQMC;

				stdErrorsMc[m][s] = stdErrMC;
				stdErrorsMcFinmath[m][s] = stdErrMcFm;
				stdErrorsQMc[m][s] = stdErrQMC;

				System.out.println("ANA: " + pricesAnalytical[m * strikes.length + s] + "\t"
						+ "MC: " + pricesMC[m][s] + " (" + stdErrMC + ")" + "\t"
						+ "MC finmath: " + pricesMCFinmath[m][s] + " (" + stdErrMcFm + ")" + "\t"
						+ "QMC: " + pricesQMC[m][s] + " (" + stdErrQMC + ")" + "\t"
						+ "BS: " + pricesBS[m][s]);
			}

		}
		/**
		 * Print prices
		 */
		System.out.println("Analytical prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesAnalytical[m * strikes.length + s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("MC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesQMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("BS prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesBS[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("Mc prices finmath");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(pricesMCFinmath[m][s] + "\t");
			}
			System.out.print("\n");
		}

		/**
		 * Print standard error
		 */
		System.out.println("MC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(stdErrorsMc[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(stdErrorsQMc[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("MC finmath std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < strikes.length; s++) {
				System.out.print(stdErrorsMcFinmath[m][s] + "\t");
			}
			System.out.print("\n");
		}
	}

	@Test
	public void getImpliedVolsForHestonTarget() throws Exception {
		/*setTargetSurfaceHeston();*/
		setTargetSurfaceBTC();

		// Get option values
		int numStrikesPerMaturity = volatilityPointsSurface.getNumberOfVolatilityPoints() / maturityGrid.length;

		double[][] pricesMC = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] pricesQMC = new double[maturityGrid.length][numStrikesPerMaturity];
		// Get analytical prices
		StopWatch sw = StopWatch.createStarted();
		double[] pricesAnalytical = lnsvqdModelAnalyticalPricer
				.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null)
				.getVolatilityPoints()
				.stream()
				.mapToDouble(volPoints -> volPoints.getVolatility())
				.toArray();
		sw.stop();
		System.out.println("time: " + sw.getTime()); // formatted string like "12.3 ms"

		double[][] stdErrorsMc = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] stdErrorsQMc = new double[maturityGrid.length][numStrikesPerMaturity];

		random.setSeed(7843657);
		for(int m = 0; m < maturityGrid.length; m++) {
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				VolatilityPoint volatilityPoint = volatilityPointsSurface.getVolatilityPoints().get(m * numStrikesPerMaturity + s);
				LocalDate date = volatilityPoint.getDate();
				double maturity = dayCountConvention.getDaycountFraction(valuationDate, date);
				double strike = volatilityPoint.getStrike();
				double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
						maturity, (int) (Math.round(maturity * 365.) * 2));


				List<Integer> seeds = random.ints(7).boxed().collect(Collectors.toList());
				double[] prices = new double[seeds.size()];
				double[] pricesQ = new double[seeds.size()];

				for(int seed : seeds) {
					// Normal MC
					// StopWatch sw1 = StopWatch.createStarted();
					LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdCallPriceSimulator.precalculatePaths(seed);
					// sw1.stop();
					// System.out.println("time MC: " + sw1.getTime());
					double simulatedOptionPrice;
					try {
						// sw1.reset();
						// sw1.start();
						simulatedOptionPrice = lnsvqdCallPriceSimulator.getImpliedVol(strike, maturity);
						// sw1.stop();
						// System.out.println("time MC: " + sw1.getTime());
					} catch(AssertionError e) {
						// System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPrice = 1000000;
					}

					prices[seeds.indexOf(seed)] = simulatedOptionPrice;

					// QMC
					LNSVQDPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
					lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
					double simulatedOptionPriceQMC;
					try {
						simulatedOptionPriceQMC = lnsvqdPriceSimulatorQMC.getImpliedVol(strike, maturity);
					} catch(AssertionError e) {
						System.err.println("Caught AssertionError: " + e.getMessage());
						simulatedOptionPriceQMC = 1000000;
					}
					pricesQ[seeds.indexOf(seed)] = simulatedOptionPriceQMC;
				}

				double averagePrice = Arrays.stream(prices).average().getAsDouble();
				double varMC = Arrays.stream(prices).map(x -> Math.pow(x - averagePrice, 2)).sum() / (seeds.size() - 1);
				double stdErrMC = Math.sqrt(varMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalMC = LNSVQDUtils.getConfidenceInterval(prices, 0.05);

				double averagePriceQMC = Arrays.stream(pricesQ).average().getAsDouble();
				double varQMC = Arrays.stream(pricesQ).map(x -> Math.pow(x - averagePriceQMC, 2)).sum() / (seeds.size() - 1);
				double stdErrQMC = Math.sqrt(varQMC) / Math.sqrt(seeds.size());
				double[] confidenceIntervalQMC = LNSVQDUtils.getConfidenceInterval(pricesQ, 0.05);

				pricesMC[m][s] = averagePrice;
				pricesQMC[m][s] = averagePriceQMC;

				stdErrorsMc[m][s] = stdErrMC;
				stdErrorsQMc[m][s] = stdErrQMC;

				System.out.println(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t"
						+ pricesMC[m][s] + "\t" + stdErrMC + "\t" + confidenceIntervalMC[0] + "\t" + confidenceIntervalMC[1] + "\t"
						+ pricesQMC[m][s] + "\t" + stdErrQMC + "\t" + confidenceIntervalQMC[0] + "\t" + confidenceIntervalQMC[1] + "\t");
			}
		}
		/**
		 * Print prices
		 */
		System.out.println("Analytical prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("MC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesQMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		/**
		 * Print standard error
		 */
		System.out.println("MC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(stdErrorsMc[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(stdErrorsQMc[m][s] + "\t");
			}
			System.out.print("\n");
		}
	}

	@Test
	public void getImpliedVolsForHestonTargetNewStdErrors() throws Exception {
		setTargetSurfaceHeston();

		// Get option values
		int numStrikesPerMaturity = volatilityPointsSurface.getNumberOfVolatilityPoints() / maturityGrid.length;

		double[][] pricesMC = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] pricesQMC = new double[maturityGrid.length][numStrikesPerMaturity];
		// Get analytical prices
		double[] pricesAnalytical = lnsvqdModelAnalyticalPricer
				.getImpliedVolSurfaceFromVolSurface(volatilityPointsSurface, null)
				.getVolatilityPoints()
				.stream()
				.mapToDouble(volPoints -> volPoints.getVolatility())
				.toArray();

		double[][] stdErrorsMc = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] stdErrorsQMc = new double[maturityGrid.length][numStrikesPerMaturity];

		for(int m = 0; m < maturityGrid.length; m++) {
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				VolatilityPoint volatilityPoint = volatilityPointsSurface.getVolatilityPoints().get(m * numStrikesPerMaturity + s);
				LocalDate date = volatilityPoint.getDate();
				double maturity = dayCountConvention.getDaycountFraction(valuationDate, date);
				double strike = volatilityPoint.getStrike();
				double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
						maturity, (int) (Math.round(maturity * 365.) * 1));

				int seed = random.nextInt();
				double implVolMC = 0;
				double implVolQMC = 0;

				double[] confidenceIntervalMC = new double[2];
				double[] confidenceIntervalQMC = new double[2];

				double stdErrMC = 0;
				double stdErrQMC = 0;

				// Normal MC
				LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
				lnsvqdCallPriceSimulator.precalculatePaths(seed);
				try {
					double[] priceStdErrorAndBounds = lnsvqdCallPriceSimulator.getPriceStdErrorAndBounds(strike, maturity);
					double price = priceStdErrorAndBounds[0];
					implVolMC = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, price);
					stdErrMC = priceStdErrorAndBounds[1];
					double lowerPrice = Math.max(price - 1.96 * stdErrMC / Math.sqrt(numberOfPaths), 1E-10);
					double upperPrice = price + 1.96 * stdErrMC / Math.sqrt(numberOfPaths);
					/*confidenceIntervalMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[2]);
					confidenceIntervalMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[3]);*/
					confidenceIntervalMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, lowerPrice);
					confidenceIntervalMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, upperPrice);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
				}

				// QMC
				int numberOfPathsQMC = 100000;
				LNSVQDPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPathsQMC, timeGrid, false);
				lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
				try {
					double[] priceStdErrorAndBounds = lnsvqdPriceSimulatorQMC.getPriceStdErrorAndBounds(strike, maturity);
					double price = priceStdErrorAndBounds[0];
					implVolQMC = lnsvqdPriceSimulatorQMC.getImpliedVolFromPrice(strike, maturity, price);
					stdErrQMC = priceStdErrorAndBounds[1];
					double lowerPrice = Math.max(price - 1.96 * stdErrQMC / Math.sqrt(numberOfPathsQMC), 1E-15);
					double upperPrice = price + 1.96 * stdErrQMC / Math.sqrt(numberOfPathsQMC);
					/*confidenceIntervalMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[2]);
					confidenceIntervalMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[3]);*/
					confidenceIntervalQMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, lowerPrice);
					confidenceIntervalQMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, upperPrice);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
				}


				System.out.println(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t"
						+ implVolMC + "\t" + stdErrMC + "\t" + confidenceIntervalMC[0] + "\t" + confidenceIntervalMC[1] + "\t"
						+ implVolQMC + "\t" + stdErrQMC + "\t" + confidenceIntervalQMC[0] + "\t" + confidenceIntervalQMC[1] + "\t");
			}
		}
		/**
		 * Print prices
		 */
		System.out.println("Analytical prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("MC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesQMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		/**
		 * Print standard error
		 */
		System.out.println("MC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(stdErrorsMc[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(stdErrorsQMc[m][s] + "\t");
			}
			System.out.print("\n");
		}
	}

	@Test
	public void testBitcoin() throws Exception {
		ArrayList<Pair<Double, Double>> strikeMatPairs = new ArrayList<>();
		/*45000.,  48000.,  55000.,  58000.,  64000.,  65000.,  70000.,
				75000.,  80000.,  85000.,  90000., 100000., 120000.]),)*/
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 45000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 48000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 55000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 58000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 64000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 65000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 70000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 75000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 80000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 85000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 90000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 100000.));
		strikeMatPairs.add(new Pair<Double, Double>(0.10122575874485597, 120000.));

		// Get option values
		int numStrikesPerMaturity = strikeMatPairs.size();

		double[][] pricesMC = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] pricesQMC = new double[maturityGrid.length][numStrikesPerMaturity];
		// Get analytical prices
		double[] pricesAnalytical = lnsvqdModelAnalyticalPricer
				.getImpliedVolsStrikeMatList(strikeMatPairs, null);

		double[][] stdErrorsMc = new double[maturityGrid.length][numStrikesPerMaturity];
		double[][] stdErrorsQMc = new double[maturityGrid.length][numStrikesPerMaturity];

		for(int m = 0; m < maturityGrid.length; m++) {
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				double maturity = strikeMatPairs.get(m * numStrikesPerMaturity + s).getKey();
				double strike = strikeMatPairs.get(m * numStrikesPerMaturity + s).getValue();
				double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
						maturity, (int) (Math.round(maturity * 365.) * 1));

				int seed = random.nextInt();
				double implVolMC = 0;
				double implVolQMC = 0;

				double[] confidenceIntervalMC = new double[2];
				double[] confidenceIntervalQMC = new double[2];

				double stdErrMC = 0;
				double stdErrQMC = 0;

				// Normal MC
				LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, true);
				lnsvqdCallPriceSimulator.precalculatePaths(seed);
				try {
					double[] priceStdErrorAndBounds = lnsvqdCallPriceSimulator.getPriceStdErrorAndBounds(strike, maturity);
					double price = priceStdErrorAndBounds[0];
					implVolMC = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, price);
					stdErrMC = priceStdErrorAndBounds[1];
					double lowerPrice = Math.max(price - 1.96 * stdErrMC / Math.sqrt(numberOfPaths), 1E-10);
					double upperPrice = price + 1.96 * stdErrMC / Math.sqrt(numberOfPaths);
					/*confidenceIntervalMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[2]);
					confidenceIntervalMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[3]);*/
					confidenceIntervalMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, lowerPrice);
					confidenceIntervalMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, upperPrice);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
				}

				// QMC
				int numberOfPathsQMC = 100000;
				LNSVQDPriceSimulatorQMC lnsvqdPriceSimulatorQMC = new LNSVQDPriceSimulatorQMC(lnsvqdModelAnalyticalPricer, numberOfPathsQMC, timeGrid, false);
				lnsvqdPriceSimulatorQMC.precalculatePaths(seed);
				try {
					double[] priceStdErrorAndBounds = lnsvqdPriceSimulatorQMC.getPriceStdErrorAndBounds(strike, maturity);
					double price = priceStdErrorAndBounds[0];
					implVolQMC = lnsvqdPriceSimulatorQMC.getImpliedVolFromPrice(strike, maturity, price);
					stdErrQMC = priceStdErrorAndBounds[1];
					double lowerPrice = Math.max(price - 1.96 * stdErrQMC / Math.sqrt(numberOfPathsQMC), 1E-15);
					double upperPrice = price + 1.96 * stdErrQMC / Math.sqrt(numberOfPathsQMC);
					/*confidenceIntervalMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[2]);
					confidenceIntervalMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, priceStdErrorAndBounds[3]);*/
					confidenceIntervalQMC[0] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, lowerPrice);
					confidenceIntervalQMC[1] = lnsvqdCallPriceSimulator.getImpliedVolFromPrice(strike, maturity, upperPrice);
				} catch(AssertionError e) {
					System.err.println("Caught AssertionError: " + e.getMessage());
				}


				System.out.println(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t"
						+ implVolMC + "\t" + stdErrMC + "\t" + confidenceIntervalMC[0] + "\t" + confidenceIntervalMC[1] + "\t"
						+ implVolQMC + "\t" + stdErrQMC + "\t" + confidenceIntervalQMC[0] + "\t" + confidenceIntervalQMC[1] + "\t");
			}
		}
		/**
		 * Print prices
		 */
		System.out.println("Analytical prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesAnalytical[m * numStrikesPerMaturity + s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("MC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC prices");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(pricesQMC[m][s] + "\t");
			}
			System.out.print("\n");
		}

		/**
		 * Print standard error
		 */
		System.out.println("MC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(stdErrorsMc[m][s] + "\t");
			}
			System.out.print("\n");
		}

		System.out.println("QMC std.errors");
		for(int m = 0; m < maturityGrid.length; m++) {
			double maturity = maturityGrid[m];
			System.out.print(maturity + "\t");
			for(int s = 0; s < numStrikesPerMaturity; s++) {
				System.out.print(stdErrorsQMc[m][s] + "\t");
			}
			System.out.print("\n");
		}
	}
	@Test
	public void printTransformedVolPathMoments() throws CalculationException {
		int numberOfPaths = 100000;
		// Get option values
		double maturity = 1.2;
		double[] timeGrid = LNSVQDUtils.createTimeGrid(0.,
				maturity, (int) Math.round(maturity * 365.));

		int seed = 1;

		LNSVQDCallPriceSimulator lnsvqdCallPriceSimulator = new LNSVQDCallPriceSimulator(lnsvqdModelAnalyticalPricer, numberOfPaths, timeGrid, false);
		lnsvqdCallPriceSimulator.precalculatePaths(seed);
		double[][][] transformedPaths = lnsvqdCallPriceSimulator.getTransformedPath();

		for(int j = 0; j < timeGrid.length; j++) {
			double moment1 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 1)).average().getAsDouble();
			double moment2 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 2)).average().getAsDouble();
			double moment3 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 3)).average().getAsDouble();
			double moment4 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 4)).average().getAsDouble();
			double moment5 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 5)).average().getAsDouble();
			double moment6 = Arrays.stream(transformedPaths[1][j]).map(x -> Math.pow(x, 6)).average().getAsDouble();
			System.out.println(timeGrid[j] +
					"\t" + moment1 +
					"\t" + moment2 +
					"\t" + moment3 +
					"\t" + moment4 +
					"\t" + moment5 +
					"\t" + moment6);
		}
	}
}