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
			0.8376,
			3.1844,
			3.058,
			1.0413,
			0.1514,
			1.8458
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

	static double[] paramInitDAXSepKappa2Positive = new double[]{
			0.137310669,
			3.1565426287720286,
			1.9367244258255474,
			0.146643909,
			-1.489064245,
			1.535139247
	};

	static double[] paramInitDAXMarch = new double[]{
			0.751963407
			, 3.560176604850655 // 3.1565426287720286,
			, 0 // 1.9367244258255474,
			, 0.467763099
			, -0.652852033
			, 0.859636078
	};

	static double[] paramInitDAXMarchKappa2Positive = new double[]{
			1
			, 4.403227439859688 // , 3.1565426287720286 // From python
			, 3.536886176176098 // , 1.5132879713105358 //, 1.9367244258255474 // From python
			, 0.2169627604050273
			, -1.4891723022148557
			, 0.9240577536198922
	};

	static double[] paramInitDAXFeb = new double[]{
			0.1598567708165421,
			3.560176604850655, // 3.1565426287720286, // From python
			0, // 1.9367244258255474, // From python
			0.189967208,
			-0.40790656321439606,
			0.20752601856648067
	};

	static double[] paramInitDAXFebKappa2Positive = new double[]{
			0.1598567708165421,
			3.1565426287720286, // From python
			1.9367244258255474, // From python
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

	static double[] paramInitDAXSepCalibKappa2Positive = new double[]{
			0.14250090891079703,
			3.1565426287720286,
			1.9367244258255474,
			0.14854371929914006,
			-1.4287502987787102,
			1.5171087093484412
	};

	static double[] paramInitDAXMarchCalib = new double[]{
			0.7183550169994655
			, 3.560176605
			, 0
			, 0.2169627604050273
			, -1.4891723022148557
			, 0.9240577536198922
	};

	static double[] paramInitDAXMarchCalibKappa2Positive = new double[]{
			0.7218799691090293,
			3.1565426287720286, // From python
			1.9367244258255474, // From python
			0.2424508550509456,
			-1.6202793135618887,
			1.151888849105794,
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

	static double[] paramInitDAXFebCalibKappa2Positive = new double[]{
			0.1666270357888853
			, 3.1565426287720286
			, 1.9367244258255474
			, 0.14281672921287228
			, -1.3894496216581893
			, 1.6057516010167643
	};

}
