#pragma once
#include <boost/multiprecision/cpp_dec_float.hpp>

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> > cpp_dec_float_100;
typedef std::vector<std::vector<cpp_dec_float_100>> matrix;
typedef std::vector<cpp_dec_float_100> vector;

int N = 13;
//cpp_dec_float_100 G(3.65);
cpp_dec_float_100 G(6.024680040776729583740234375);
//cpp_dec_float_100 G(5.15);
//cpp_dec_float_100 e(2.7182818284590452353602874713526624977572470937);


/* calculate combination number */
cpp_dec_float_100 combination(cpp_dec_float_100 n, cpp_dec_float_100 m)
{
	/* Cnm = Cn(n - m)*/
	m = m > n - m ? n - m : m;
	cpp_dec_float_100 result = 1;
	cpp_dec_float_100 loops = 1;
	while (loops <= m) {
		result *= n--;
		result /= loops++;
	}
	return result;
}

cpp_dec_float_100 factorial(cpp_dec_float_100 n)
{
	cpp_dec_float_100 result = 1.0;
	while (n)
		result *= n--;
	return result;
}

/* generate matrix B, n is the dimension*/
void generateB(matrix & b, int n)
{
	std::cout << "B matrix" << std::endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == 0)
				b[i][j] = 1;
			else if (i > 0 && j >= i)
				b[i][j] = std::pow(-1, j - i) * combination(i + j - 1, j - i);
			else
				b[i][j] = 0;
			std::cout << "\t" << b[i][j];
		}
		std::cout << std::endl;
	}
}

void generateC(matrix & c, int n)
{
	std::cout << "C matrix" << std::endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i < j)
				c[i][j] = 0;
			else if (i == j && i == 0)
				c[i][j] = 0.5;
			else {
				//c[i][j] = std::pow(-1, i - j) * boost::multiprecision::pow((cpp_dec_float_100)4, j) * i * 
				//	factorial(i + j - 1) / (factorial(i - j) * factorial(2 * j));
				for (int k = 0; k <= i; ++k)
					c[i][j] += combination(2 * i, 2 * k) * combination(k, k + j - i);
				c[i][j] *= std::pow(-1, i - j);
			}
			std::cout << "\t" << c[i][j];
		}
		std::cout << std::endl;
	}
}

void generateDr(matrix & dr, int n)
{
	std::cout << "D matrix" << std::endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i != j)
				dr[i][j] = 0;
			else if (i == j && i == 0)
				dr[i][j] = 1;
			else if (i == j && j == 1)
				dr[i][j] = -1;
			else
				dr[i][j] = dr[i - 1][j - 1] * 2 * (2 * i - 1) / (i - 1);
			std::cout << "\t" << dr[i][j];
		}
		std::cout << std::endl;
	}
}

void generateF(vector & f, int n)
{
	std::cout << "F vector" << std::endl;
	for (int i = 0; i < n; ++i) {
		f[i] = factorial(2 * i) * boost::multiprecision::exp((cpp_dec_float_100)i + 0.5L) 
			/ (factorial(i) * boost::multiprecision::pow((cpp_dec_float_100)2, 2 * i - 1) 
				* boost::multiprecision::pow((cpp_dec_float_100)i + G + 0.5L, (cpp_dec_float_100)i + 0.5L));
		std::cout << "\t" << f[i];
	}
	std::cout << std::endl;

}

//template<MatrixT>
matrix matrixMultiply(matrix & a, matrix & b)
{
	matrix result(a.size(), vector (a[0].size()));
	for (int i = 0; i < a.size(); ++i) {
		for (int j = 0; j < b[0].size(); ++j) {
			for (int k = 0; k < b[0].size(); ++k)
				result[i][j] += a[i][k] * b[k][j];
			//std::cout << "\t" << result[i][j];
		}
		//std::cout << std::endl;
	}
	return result;
}

vector matrixVectorMultiply(matrix & a, vector & b)
{
	std::cout << "result" << std::endl;
	vector result(b.size());
	for (int i = 0; i < a.size(); ++i) {
		for (int k = 0; k < b.size(); ++k)
			result[i] += a[i][k] * b[k];
		//std::cout << std::setprecision(10) << "\t" << result[i];
	}
	
	return result;
}

vector generateA()
{
	int n = N;
	matrix B(n, vector(n));
	generateB(B, n);
	matrix C(n, vector(n));
	generateC(C, n);
	matrix Dr(n, vector(n));
	generateDr(Dr, n);
	vector F(n);
	generateF(F, n);
	matrix r = matrixMultiply(Dr, B);
	r = matrixMultiply(r, C);
	vector a = matrixVectorMultiply(r, F);
	for (int i = 0; i < n; ++i)
		std::cout << std::setprecision(10) << "\t" << a[i];
	//cpp_dec_float_100 pi = boost::math::constants::pi<cpp_dec_float_100>();
	//cpp_dec_float_100 W(boost::multiprecision::exp(G) / boost::multiprecision::sqrt(2 * pi));

	//vector c(n);
	//for (int i = 0; i < n; ++i) {
	//	c[i] = a[i] * W;
	//	std::cout << std::setprecision(10) << "\t" << c[i];
	//}
	//return c;
	return a;
}

template <class Accur = double>
inline Accur lanzocsConstWithE(Accur z)
{
	double a[] = {
		0.00606184234624890652578375396455593688322246366549701353315381,
		1.42562835481480388911205082410391681869770004663319327370524,
		-2.14753419549590132870001860978495612225150976745591920430031,
		0.957266911876239916727464187654754930799056796426414254410443,
		-0.128688659284825666157814588986354869299565554291793485197464,
		0.00308864347836427206444800003150217853189566711517579490570701,
		-9.7849610972063267859543785146578805638938457810257219936838e-07,
		-1.77688484056248341590450071107552867435518550688147157442928e-08,
		1.98515619784786124336217205348896936351743198446648459440374e-08,
		-1.24774538787442127470155013366827607159775465966278217150944e-08,
		5.60962322320093594808388820009754037031162213257143056660105e-09,
		-1.61332427088018334077343069386684131191662695730605030781635e-09,
		2.19109779267020750367096552334711876334150700512296792774128e-10
	};

	double s = a[0];
	for (int i = 1; i < N; ++i)
		s += a[i] / (z + i);

	return s;
}

//template <class Accur = double>
//inline Accur lanzocsConstWithE2(Accur z)
//{
//	int n = N;
//	double result = a[0];
//	double zserial = 1.0;
//	for (int i = 1; i < n; ++i) {
//		zserial *= z + i;
//	}
//
//	double zasum = 0;
//	for (int i = 1; i < n; ++i) {
//		zasum += a[i] * zserial / (z + i);
//	}
//	result = (result * zserial + zasum) / zserial;
//	return result;
//}

