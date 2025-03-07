package net.finmath.equities;

import net.finmath.functions.AnalyticFormulas;
import net.finmath.time.daycount.DayCountConvention_ACT_365;
import org.junit.Assert;
import org.junit.Test;

import java.time.LocalDate;

import static java.lang.Math.E;

public class GeneralTests {

	@Test
	public void comparaBSPriceAndImppliedVol() {
		LocalDate valuationDate = LocalDate.parse("2024-09-30");
		LocalDate maturity = LocalDate.parse("2025-06-20");
		DayCountConvention_ACT_365 dCC = new DayCountConvention_ACT_365();

		double spot = 1;
		double ttm = dCC.getDaycountFraction(valuationDate, maturity);
		double r = 0;
		double forward = 1;
		double strike = 1.4;

		double vol1 = 0.1256211050901605;
		double vol2 = 0.5;

		double bsPrice1 = 5.8518458685649172E-5; // AnalyticFormulas.blackScholesOptionValue(spot, r, vol1, ttm, strike);
		double bsPrice2 = AnalyticFormulas.blackScholesOptionValue(spot, r, vol2, ttm, strike);

		double impliedVol1 = AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, ttm, strike, 1, bsPrice1);
		double impliedVol2 = AnalyticFormulas.blackScholesOptionImpliedVolatility(forward, ttm, strike, 1, bsPrice2);

		Assert.assertEquals(vol1 , impliedVol1, E-3);
		// Assert.assertEquals(vol2 , impliedVol2, E-3);
	}

	@Test
	public void test() {
		double[] arr = new double[]{1, 2, 3};
		double p = arr[1];
		System.out.println(p);
		arr[1] = 0;
		System.out.println(p);
	}
}
