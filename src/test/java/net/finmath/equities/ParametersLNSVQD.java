package net.finmath.equities;

public class ParametersLNSVQD {
	static double[] paramVectorHestonDAXSep = new double[]{
			0.024536987
			, 4.
			, 0.036676304
			, 1.211288333
			, -0.672024524
	};

	static double[] paramVectorHestonDAXMarch = new double[]{
			0.524442119
			, 4
			, 0.086423244
			, 2.190580015
			, -0.73124787
	};

	static double[] paramVectorHestonDAXFeb = new double[]{
			0.034583515
			, 4
			, 0.03608754
			, 1.2835035
			, -0.634272101
	};
	static double[] paramVectorBitcoin = new double[]{
			5, // 0.8376,
			4.8606,
			4.7938,
			1.0139,
			0.1985,
			2.3690
	};

	static double[] paramBlackScholes = new double[]{
			0.2,
			0,
			0,
			0,
			0,
			0
	};

	static double[] paramInitDAXSep = new double[]{
			0.137310669,
			3.560176604850655, // 3.1565426287720286,
			0, // 1.9367244258255474,
			0.146643909,
			-1.489064245,
			1.535139247
	};

	static double[] paramInitDAXMarch = new double[]{
			0.751963407
			, 3.1565426287720286 // From python
			, 1.9367244258255474 // From python
			, 0.467763099
			, -0.652852033
			, 0.859636078
	};

	static double[] paramInitDAXFeb = new double[]{
			0.1598567708165421,
			3.560176604850655, // 3.1565426287720286, // From python
			0, // 1.9367244258255474, // From python
			0.189967208,
			-0.40790656321439606,
			0.20752601856648067
	};

	static double[] paramInitDAXSepCalib = new double[]{
			0.14073129923299094,
			3.560176604850655,
			0,
			0.1459432399874267,
			-1.41428038221871,
			1.4653117049022548
	};

	static double[] paramInitDAXMarchCalibOld = new double[]{
			0.78096181461525,
			3.1565426287720286, // From python
			1.9367244258255474, // From python
			0.4961670369028076,
			-0.7351529740014809,
			0.6002284343654938,
	};

	static double[] paramInitDAXMarchCalib = new double[]{
			0.7183550169994655
			, 3.560176605
			, 0
			, 0.2169627604050273
			, -1.4891723022148557
			, 0.9240577536198922
	};

	static double[] paramInitDAXFebCalibOld = new double[]{
			0.166627036,
			3.156542629,
			1.936724426,
			0.142816729,
			-1.389449622,
			1.605751601,
	};

	static double[] paramInitDAXFebCalib = new double[]{
			0.162735884
			, 3.560176605
			, 0
			, 0.142510755
			, -1.372319642
			, 1.54799513
	};

}
