/*
言語:
C++
OS:
Microsoft Windows 10 Home (64 bit)
コンパイル:
C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall x86_amd64 & cl /EHsc /O2 /nologo $file_name
g++ -O2 --exec-charset=CP932 -std=c++14 $file_name -o $file_base_name
日付:
2017/01/12
*/
#ifndef MYUTIL_H
#define MYUTIL_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <sstream>		// ss
#include <fstream>		// ifstream, ofstream
#include <iomanip>		// マニピュレータ
#include <string>
#include <vector>
#include <time.h>

// Eigen 行列計算
#define EIGEN_NO_DEBUG
#include <Eigen/Core>

// ########################################
// マクロ
// Repetition
#define FOR(i,a,b)		for(int i=(a);i<(b);++i)
#define REP(i,n)		FOR(i,0,n)
#define RFOR(i,a,b)		for(int i=(a);i>=(b);--i)

namespace MyUtil {

using namespace std;
const std::ios::fmtflags IOMANIP_FLAGS_SAVED = std::cout.flags();
const double PI = M_PI;
const int MY_INT_MAX = 2147483647;
const int MY_INT_MIN = -2147483647-1;
using status_t = int;		// 成功失敗などの返り値用
const status_t STATUS_SUCCESS = 1;
const status_t STATUS_FAILURE = 0;

// 汎用Print
template <typename T> void p(T a);
template <typename T> void pn(T a);
void _ps();
template <typename T, typename... ARGV> void _ps(T first, ARGV... argv);
template <typename... ARGV> void ps(ARGV... argv);
template <typename... ARGV> void psn(ARGV... argv);
// void _psd(const string delimiter);
// template <typename T, typename... ARGV> void _psd(const string delimiter, T first, ARGV... argv);
template <typename... ARGV> void psd(const string delimiter, ARGV... argv);
template <typename... ARGV> void psdn(const string delimiter, ARGV... argv);
// template <typename T> void pa(T a[], int size);
// template <typename T> void pan(T a[], int size);
template <typename T> void PrintV(const vector<T>& v, const string delimiter = "\t");
template <typename T> void PrintVV(const vector<vector<T>>& vv, const string delimiter = "\t");
// 汎用数学関数
double Rad2Deg(double rad);
double Deg2Rad(double deg);
double NormalizeRad(double input);
double NormalizeRadSym(double input);
template<class T> inline T sq(T x);
// template <typename T> void SwapVector(vector<T>& v, int row1, int row2);
// template <typename T> void SwapVector2d(vector<vector<T>>& vv, int col1, int col2);
// template <typename T> void SwapVector2dRow(vector<vector<T>>& vv, int row1, int row2);
// template <typename T> void SwapVector2dCol(vector<vector<T>>& vv, int col1, int col2);
// アルゴリズム
// template <typename T1, typename T2> T1 Nibutan(T2 (*f)(T1 x), T1 min, T1 max, double epsilon);
// double Nibutan(double (*f)(double x, double t), double min, double max, double epsilon, double t);
// 数値計算
// void RungeKutta(double (*f[])(double t, double *x), double *x, double t0, double tn, double dt, int num);
// double GaussJordanWithPivoting(vector<vector<double>>& A, vector<double>& b, int size);
// int _GaussJordanGetMaxRow(vector<vector<double>>& A, int size, int minRow);
// void _GaussJordanSwapRow(vector<vector<double>>& A, vector<double>& b, int size, int row1, int row2);
// int Jacobi(vector<vector<double>>& A, vector<double>& b, int size, double epsilon, string outputFileName);
// void _JacobiOutputStdout(int count, Eigen::VectorXd& ex, double residualNorm, int size);
// void _JacobiOutputFile(int count, Eigen::VectorXd& ex, double residualNorm, int size, fstream& fs);
// 行列計算Eigen
void Vector2d2EigenMatrix2d(vector<vector<double>>& A, Eigen::MatrixXd& eA);
void Vector1d2EigenVector1d(vector<double>& a, Eigen::VectorXd& ea);
int EigenMatrix2d2Vector2d(vector<vector<double>>& A, Eigen::MatrixXd& eA);
int EigenVector1d2Vector1d(vector<double>& a, Eigen::VectorXd& ea);
// int TransposeEigenMatrix2d(Eigen::MatrixXd& eA);
// double GetMaxEigenVector1dError(Eigen::VectorXd& epreX, Eigen::VectorXd& ex, int size);
// 文字列処理系
string ReplaceString(const string& target, const string& from, const string& to);
vector<string> SplitString(const string& input, const string delimiter);
vector<string> SplitStringSs(const string& input, const string delimiter=" ");
template <typename T>
vector<T> CastStringVector(const vector<string>& input);
template <typename T>
string ToString(const T& input);
// Util
string GetStringTime();
inline status_t outputFile(const string& filename, const string& str, int isApp = 0);
inline status_t outputFile(const string& filename, const ostringstream& oss, int isApp = 0);
template <typename T>
status_t output2dVector2Csv(const string& filename, vector<vector<T>>& vv);
template <typename T>
status_t output1dVector2Csv(const string& filename, vector<T>& v, int isRow = 1);
template <typename T>
string Matrix2Gnuplot3DGridData(vector<vector<T>>& vv, double dx, double dy, int stepX, int stepY, double offsetX=0, double offsetY=0);
template <typename T>
string Matrix2GnuplotPm3dMapData(vector<vector<T>>& vv, double dx, double dy, double offsetX=0, double offsetY=0);
// string GetBar(int num);
string GetDbar(int num);
string GetENotation(double num, int precision, int isPos=0, int exponentDigit=2);
string GetFNotation(double num, int precision, int isPos=0);
ios::fmtflags SetIomanipScientificPos(int precision);
void ResetIomanip();
void _PrintErr(const string& s);


// ########################################
// 汎用Print
template <typename T>
void p(T a) {
	cout << a << flush;
}

template <typename T>
void pn(T a) {
	cout << a << endl;
}

void _ps() {
}

template <typename T, typename... ARGV>
void _ps(T first, ARGV... argv) {
	// cout << first << " ";
	cout << first;
	_ps(argv...);
}

template <typename... ARGV>
void ps(ARGV... argv) {
	_ps(argv...);
	cout << flush;
}

template <typename... ARGV>
void psn(ARGV... argv) {
	_ps(argv...);
	cout << endl;
}

void _psd(const string delimiter) {
}

template <typename T, typename... ARGV>
void _psd(const string delimiter, T first, ARGV... argv) {
	cout << first << delimiter;
	_psd(delimiter, argv...);
}

template <typename... ARGV>
void psd(const string delimiter, ARGV... argv) {
	_psd(delimiter, argv...);
	cout << flush;
}

template <typename... ARGV>
void psdn(const string delimiter, ARGV... argv) {
	_psd(delimiter, argv...);
	cout << endl;
}

// template <typename T>
// void pa(T a[], int size) {
// 	for (int i=0; i<size; ++i) {
// 		cout << a[i] << " ";
// 	}
// }

// template <typename T>
// void pan(T a[], int size) {
// 	pa<T>(a, size);
// 	cout << endl;
// }

template <typename T>
void PrintV(const vector<T>& v, const string delimiter) {
	for (auto v2: v) {
		cout << v2 << delimiter;
	}
	cout << endl;
}

template <typename T>
void PrintVV(const vector<vector<T>>& vv, const string delimiter) {
	for (auto v1: vv) {
		for (auto v2: v1) {
			cout << v2 << delimiter;
		}
		cout << endl;
	}
}


// ########################################
// 汎用数学関数
double Rad2Deg(double rad) {
	return rad * 180 / PI;
}

double Deg2Rad(double deg) {
	return deg * PI / 180.0;
}

//radianを0~2PIにする
double NormalizeRad(double input) {
	if (input >= 0 ) {
		return input - (int)(input / (2*PI)) *2*PI;
	} else {
		return input + ((int)((-1*input) / (2*PI)) + 1)* 2*PI;
	}
}

//radianを-PI~PIにする
double NormalizeRadSym(double input) {
	double result = NormalizeRad(input);
	if (result > PI) {
		result = result - 2*PI;
	}
	return result;
}

template<class T>
inline T sq(T x) {return x*x;}

// template <typename T>
// void SwapVector(vector<T>& v, int row1, int row2) {
// 	T temp = v[row1];
// 	v[row1] = v[row2];
// 	v[row2] = temp;
// }

// template <typename T>
// void SwapVector2d(vector<vector<T>>& vv, int col1, int col2) {
// 	for (int i=0; i<vv.size(); ++i) {
// 		SwapVector<T>(vv[i], col1, col2);
// 	}
// }

// template <typename T>
// void SwapVector2dRow(vector<vector<T>>& vv, int row1, int row2) {
// 	SwapVector<vector<T>>(vv, row1, row2);
// }

// template <typename T>
// void SwapVector2dCol(vector<vector<T>>& vv, int col1, int col2) {
// 	SwapVector2d<T>(vv, col1, col2);
// }


// ########################################
// アルゴリズム
// template <typename T1, typename T2>
// T1 Nibutan(T2 (*f)(T1 x), T1 min, T1 max, double epsilon) {
// 	T1 result = 0;
// 	T2 temp = 0;
// 	while ((max - min) > epsilon) {
// 		result = (max + min) / 2;
// 		temp = (*f)(result);
// 		if (temp < 0) {
// 			min = result;
// 		} else {
// 			max = result;
// 		}
// 	}
// 	return result;
// }

// double Nibutan(double (*f)(double x, double t), double min, double max, double epsilon, double t) {
// 	double result = 0;
// 	double temp = 0;
// 	while ((max - min) > epsilon) {
// 		result = (max + min) / 2.0;
// 		temp = (*f)(result, t);
// 		if (temp < 0) {
// 			min = result;
// 		} else {
// 			max = result;
// 		}
// 	}
// 	return result;
// }


// ########################################
// 数値計算
// void RungeKutta(double (*f[])(double t, double *x), double *x, double t0, double tn, double dt, int num) {
// 	vector<double> k1(num), k2(num), k3(num), k4(num);		// RungeKuttaの係数
// 	vector<double> temp(num);

// 	double t = t0;

// 	while (t < tn) {
// 		for (int j=0; j<num; ++j) {
// 			k1[j] = (*f[j])(t, x);
// 			temp[j] = x[j] + dt*k1[j]/2;
// 		}
// 		for (int j=0; j<num; ++j) {
// 			k2[j] = (*f[j])(t+dt/2, &temp.front());
// 		}
// 		for (int j=0; j<num; ++j) {
// 			temp[j] = x[j] + dt*k2[j]/2;
// 		}
// 		for (int j=0; j<num; ++j) {
// 			k3[j] = (*f[j])(t+dt/2, &temp.front());
// 		}
// 		for (int j=0; j<num; ++j) {
// 			temp[j] = x[j] + dt*k3[j];
// 		}
// 		for (int j=0; j<num; ++j) {
// 			k4[j] = (*f[j])(t+dt, &temp.front());
// 			x[j] += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])*dt/6;
// 		}
// 		t += dt;
// 	}
// }

// // A x = b を解く
// // size : 行列，ベクトルサイズ
// // 返り値はesidual 2-Norm
// double GaussJordanWithPivoting(vector<vector<double>>& A, vector<double>& b, int size) {
// 	double residualNorm;			// Residual 2-Norm
// 	// Residual 2-Norm計算のため，初期値を保存
// 	vector<vector<double>> AInit = A;
// 	vector<double> bInit = b;
// 	for (int i=0; i<size; ++i) {
// 		// Pivoting
// 		if (i < (size-1)) {
// 			int maxRow = _GaussJordanGetMaxRow(A, size, i);
// 			// cout << maxRow << "\n";
// 			_GaussJordanSwapRow(A, b, size, i, maxRow);
// 		}

// 		// PrintVV(A);
// 		// MyUtil::PrintV(b);
// 		// cout << "\n";
// 		// getchar();

// 		// 行の正規化
// 		double Aii = A[i][i];
// 		double AiiInverse = 1.0 / Aii;		// 除算回数を減らすため
// 		// 下三角成分は0のはずだが，エラーチェックのためにそこも操作する．
// 		for (int j=0; j<size; ++j) {
// 			A[i][j] *= AiiInverse;
// 		}
// 		b[i] *= AiiInverse;

// 		// PrintVV(A);
// 		// MyUtil::PrintV(b);
// 		// cout << "\n";
// 		// getchar();

// 		// 消去
// 		for (int j=0; j<size; ++j) {
// 			if (j != i) {
// 				double Aji = A[j][i];
// 				// 下三角成分は0のはずだが，エラーチェックのためにそこも操作する．
// 				for (int k=0; k<size; ++k) {
// 					A[j][k] = A[j][k] - Aji * A[i][k];
// 				}
// 				b[j] = b[j] - Aji * b[i];
// 			}
// 		}

// 		// PrintVV(A);
// 		// MyUtil::PrintV(b);
// 		// cout << "\n";
// 		// getchar();
// 	}

// 	// Residual 2-Normの計算
// 	// STLからEigenへ変換
// 	Eigen::MatrixXd eAInit;
// 	Eigen::VectorXd ebInit;
// 	Eigen::VectorXd ex;
// 	MyUtil::Vector2d2EigenMatrix2d(AInit, eAInit);
// 	MyUtil::Vector1d2EigenVector1d(bInit, ebInit);
// 	MyUtil::Vector1d2EigenVector1d(b, ex);
// 	// cout << (eAInit) << endl;
// 	// cout << (ex) << endl;
// 	// cout << (ebInit) << endl;
// 	// cout << (eAInit * ex - ebInit) << endl;
// 	residualNorm = (eAInit * ex - ebInit).norm();
// 	return residualNorm;
// }

// int _GaussJordanGetMaxRow(vector<vector<double>>& A, int size, int minRow) {
// 	int maxRow = minRow;
// 	int k = minRow;
// 	double nowMax = abs(A[k][k]);
// 	for (int i=minRow+1; i<size; ++i) {
// 		if (nowMax < abs(A[i][k])) {
// 			nowMax = A[i][k];
// 			maxRow = i;
// 		}
// 	}
// 	return maxRow;
// }

// void _GaussJordanSwapRow(vector<vector<double>>& A, vector<double>& b, int size, int row1, int row2) {
// 	SwapVector2dRow<double>(A, row1, row2);
// 	SwapVector<double>(b, row1, row2);

// /* 旧版
// 	vector<double> vtmp = A[col1];
// 	A[col1] = A[col2];
// 	A[col2] = vtmp;
// 	double dtmp = b[col1];
// 	b[col1] = b[col2];
// 	b[col2] = dtmp;
// */
// }

// int Jacobi(vector<vector<double>>& A, vector<double>& b, int size, double epsilon, string outputFileName) {
// 	using StlMatrix = vector<vector<double>>;
// 	StlMatrix L = StlMatrix(size, vector<double>(size, 0));	// 下三角
// 	StlMatrix D = StlMatrix(size, vector<double>(size, 0));	// 対角
// 	StlMatrix U = StlMatrix(size, vector<double>(size, 0));	// 上三角
// 	for (int i=0; i<size; ++i) {
// 		D[i][i] = A[i][i];
// 	}
// 	for (int i=0; i<size; ++i) {
// 		for (int j=i+1; j<size; ++j) {
// 			U[i][j] = A[i][j];
// 		}
// 	}
// 	for (int i=0; i<size; ++i) {
// 		for (int j=0; j<i; ++j) {
// 			L[i][j] = A[i][j];
// 		}
// 	}

// 	// STLからEigenへ変換
// 	Eigen::MatrixXd eA, eL, eD, eU;
// 	Eigen::VectorXd eb;
// 	Vector2d2EigenMatrix2d(A, eA);
// 	Vector2d2EigenMatrix2d(L, eL);
// 	Vector2d2EigenMatrix2d(D, eD);
// 	Vector2d2EigenMatrix2d(U, eU);
// 	Vector1d2EigenVector1d(b, eb);

// 	fstream fs;
// 	fs.open(outputFileName, ios::out);
// 	if(! fs.is_open()) {
// 		return EXIT_FAILURE;
// 	}

// 	int count = 0;
// 	const int mod = 250;
// 	// 初期値生成
// 	Eigen::VectorXd ex = Eigen::VectorXd::Zero(size);
// 	for (int i=0; i<size; ++i) {
// 		ex(i) = eb(i) / eD(i,i);
// 	}
// 	Eigen::VectorXd epreX = ex;
// 	double residualNorm = (eA * ex - eb).norm();

// 	MyUtil::p("LOOP\t");
// 	fs << "!LOOP\t";
// 	for (int i=0; i<size; ++i) {
// 		MyUtil::ps("x", i+1, "\t");
// 		fs << "x" << i+1 << "\t";
// 	}
// 	MyUtil::p("err\n");
// 	fs << "err" << "\n";
// 	MyUtil::p("================================================================\n");
// 	fs << "!================================================================" << "\n";

// 	while (1) {
// 		if (count % mod == 0) {
// 			_JacobiOutputStdout(count, ex, residualNorm, size);
// 		}
// 		_JacobiOutputFile(count, ex, residualNorm, size, fs);
// 		epreX = ex;
// 		ex = (-eL - eU) * epreX + eb;
// 		for (int i=0; i<size; ++i) {
// 			ex(i) = ex(i) / eD(i,i);
// 		}
// 		count++;
// 		// 旧版
// 		// residualNorm = GetMaxEigenVector1dError(epreX, ex, size);
// 		residualNorm = (eA * ex - eb).norm();
// 		if (residualNorm < epsilon) {
// 			break;
// 		}
// 	}

// 	MyUtil::p("================================================================\n");
// 	_JacobiOutputStdout(count, ex, residualNorm, size);
// 	_JacobiOutputFile(count, ex, residualNorm, size, fs);

// 	fs.close();
// 	std::cout.flags(IOMANIP_FLAGS_SAVED);
// 	return 1;
// }

// void _JacobiOutputStdout(int count, Eigen::VectorXd& ex, double residualNorm, int size) {
// 	MyUtil::ps("[", count, "]");
// 	for (int i=0; i<size; ++i) {
// 		int precision = ex(i) < 0 ? 4 : 5;
// 		cout << "\t" << fixed << setprecision(precision) << ex(i);
// 	}
// 	cout << "\t" << scientific << setprecision(2) << residualNorm << "\n";
// }

// void _JacobiOutputFile(int count, Eigen::VectorXd& ex, double residualNorm, int size, fstream& fs) {
// 	fs << count;
// 	for (int i=0; i<size; ++i) {
// 		fs << "\t" << fixed << setprecision(12) << ex(i);
// 	}
// 	fs << "\t" << scientific << setprecision(2) << residualNorm << "\n";
// }


// ########################################
// 行列計算Eigen
void Vector2d2EigenMatrix2d(vector<vector<double>>& A, Eigen::MatrixXd& eA) {
	int rows = A.size();
	int cols = A[0].size();
	vector<double> A1d = vector<double>(rows*cols, 0);
	for (int i=0; i<cols; ++i) {
		for (int j=0; j<rows; ++j) {
			A1d[i*rows+j] = A[j][i];
		}
	}
	eA = Eigen::Map<Eigen::MatrixXd>(&A1d.front(), rows, cols);
	// TransposeEigenMatrix2d(eA);
}

void Vector1d2EigenVector1d(vector<double>& a, Eigen::VectorXd& ea) {
	int size = a.size();
	ea = Eigen::Map<Eigen::VectorXd>(&a.front(), size);
}

int EigenMatrix2d2Vector2d(vector<vector<double>>& A, Eigen::MatrixXd& eA) {
	int rows = eA.rows();
	int cols = eA.cols();
	if (A.size() != rows) {
		_PrintErr("Not Match Rows at EigenMatrix2d2Vector2d.");
		return 0;
	}
	if (A[0].size() != cols) {
		_PrintErr("Not Match Cols at EigenMatrix2d2Vector2d.");
		return 0;
	}
	for (int i=0; i<rows; ++i) {
		for (int j=0; j<cols; ++j) {
			A[i][j] = eA(i,j);
		}
	}
	return 1;
}

int EigenVector1d2Vector1d(vector<double>& a, Eigen::VectorXd& ea) {
	int size = ea.size();
	if (a.size() != size) {
		_PrintErr("Not Match Size at EigenVector1d2Vector1d.");
		return 0;
	}
	for (int i=0; i<size; ++i) {
		a[i] = ea(i);
	}
	return 1;
}

// // 参照でメモリを消費せずに転置
// int TransposeEigenMatrix2d(Eigen::MatrixXd& eA) {
// 	int size = eA.cols();
// 	if (size != eA.rows()) {
// 		_PrintErr("Not Square Matrix at TransposeEigenMatrix2d.");
// 		return 0;
// 	}
// 	double temp;
// 	for (int i=0; i<size; ++i) {
// 		for (int j=0; j<i; ++j) {
// 			temp = eA(i,j);
// 			eA(i,j) = eA(j,i);
// 			eA(j,i) = temp;
// 		}
// 	}
// 	return 1;
// }

// double GetMaxEigenVector1dError(Eigen::VectorXd& epreX, Eigen::VectorXd& ex, int size) {
// 	double maxErr = 0;
// 	for (int i=0; i<size; ++i) {
// 		double err = abs(epreX(i) - ex(i));
// 		if (maxErr < err) {
// 			maxErr = err;
// 		}
// 	}
// 	return maxErr;
// }


// ########################################
// 文字列処理系
// string ReplaceString(const string target, const string from, const string to) {
string ReplaceString(const string& target, const string& from, const string& to) {
	string result = target;
	string::size_type pos = 0;
	while (pos = result.find(from, pos), pos != string::npos) {
		result.replace(pos, from.length(), to);
		pos += to.length();
	}
	return result;
}

vector<string> SplitString(const string& input, const string delimiter) {
	vector<string> result;
	if (input.size() == 0) return result;
	// delimiter = "" の場合
	if (delimiter == "" || delimiter.size() == 0) {
		for (int i=0;i<input.size();++i) {
			result.push_back(input.substr(i,1));
		}
		return result;
	}
	string::size_type pos = 0;
	while(pos != string::npos) {
		string::size_type p = input.find(delimiter, pos);
		if(p == string::npos) {
			result.push_back(input.substr(pos));
			break;
		} else {
			result.push_back(input.substr(pos, p - pos));
		}
		pos = p + delimiter.size();
	}
	return result;
}

vector<string> SplitStringSs(const string& input, const string delimiter) {
	string str = ReplaceString(input, delimiter, " ");
	vector<string> output;
	string elem;
	istringstream iss(str);
	while (iss >> elem) {
		output.push_back(elem);
	}
	return output;
}

template <typename T>
vector<T> CastStringVector(const vector<string>& input) {
	int size = input.size();
	vector<T> result(size);
	for (int i=0;i<size;++i) {
		istringstream iss(input[i]);
		iss >> result[i];
	}
	return result;
}

// MinGWのto_stringにはバグがあるので代替関数
template <typename T>
string ToString(const T& input) {
	ostringstream oss;
	oss << input;
	return oss.str();
}


// ########################################
// Util
string GetStringTime() {
	stringstream ss;
	int y, m, d, hor, min, sec;
	time_t now = time(NULL);
	struct tm *pnow = localtime(&now);
	// char week[][3] = {"日","月","火","水","木","金","土"};
	// printf("今日は%2d年%2d月%2d日(%s)です。\n",
	// pnow->tm_year+1900,
	// 	pnow->tm_mon + 1,
	// 	pnow->tm_mday,
	// 	week[pnow->tm_wday]);

	y = pnow->tm_year+1900;
	m = pnow->tm_mon + 1;
	d = pnow->tm_mday;
	hor = pnow->tm_hour;
	min = pnow->tm_min;
	sec = pnow->tm_sec;

	ss << setfill('0') << setw(4) << right << y << "-";
	ss << setfill('0') << setw(2) << right << m << "-";
	ss << setfill('0') << setw(2) << right << d << "_";
	ss << setfill('0') << setw(2) << right << hor << "-";
	ss << setfill('0') << setw(2) << right << min << "-";
	ss << setfill('0') << setw(2) << right << sec;
	return ss.str();
}

inline status_t outputFile(const string& filename, const string& str, int isApp) {
	ofstream ofs;
	// ofs.open(filename, ios::out);
	if (isApp == 1) {
		ofs.open(filename, ios::app);
	} else {
		ofs.open(filename);
	}
	if(! ofs.is_open()) {
		return STATUS_FAILURE;
	}
	ofs << str;
	ofs.close();
	return STATUS_SUCCESS;
}

inline status_t outputFile(const string& filename, const ostringstream& oss, int isApp) {
	return outputFile(filename, oss.str(), isApp);
}

template <typename T>
status_t output2dVector2Csv(const string& filename, vector<vector<T>>& vv) {
	ostringstream oss;
	oss.precision(14);
	int rows = vv.size();
	int cols = vv[0].size();
	for (int i=0; i<rows; ++i) {
		for (int j=0; j<cols-1; ++j) {
			oss << vv[i][j] << ",";
		}
		oss << vv[i][cols-1] << "\n";
	}
	return outputFile(filename, oss);
}

template <typename T>
status_t output1dVector2Csv(const string& filename, vector<T>& v, int isRow) {
	ostringstream oss;
	int size = v.size();
	if (isRow == 1) {
		for (int i=0; i<size-1; ++i) {
			oss << v[i] << ",";
		}
		oss << v[size-1] << "\n";
	} else {
		for (int i=0; i<size; ++i) {
			oss << v[i] << "\n";
		}
	}
	return outputFile(filename, oss);
}

template <typename T>
string Matrix2Gnuplot3DGridData(vector<vector<T>>& vv, double dx, double dy, int stepX, int stepY, double offsetX, double offsetY) {
	stringstream ss;
	int NX = vv.size();
	int NY = vv[0].size();
	for (int i=0; i<NX; i+=stepX) {
		for (int j=0; j<NY; ++j) {
			ss << i*dx+offsetX << "\t" << j*dy+offsetY << "\t" << vv[i][j] << "\n";
		}
		ss << "\n\n";
	}
	for (int i=0; i<NY; i+=stepY) {
		for (int j=0; j<NX; ++j) {
			ss << j*dx+offsetX << "\t" << i*dy+offsetY << "\t" << vv[j][i] << "\n";
		}
		ss << "\n\n";
	}
	return ss.str();
}

template <typename T>
string Matrix2GnuplotPm3dMapData(vector<vector<T>>& vv, double dx, double dy, double offsetX, double offsetY) {
	string result = Matrix2Gnuplot3DGridData(vv, dx, dy, 1, 1, offsetX, offsetY);
	return ReplaceString(result, "\n\n\n", "\n\n");
}

// string GetBar(int num) {
// 	string bar(num, '-');
// 	return bar;
// }

string GetDbar(int num) {
	string dbar(num, '=');
	return dbar;
}

string GetENotation(double num, int precision, int isPos, int exponentDigit) {
	stringstream ss;
	int exponent = (num==0) ? 0 : static_cast<int>(log10(abs(num)));
	if (exponent < 0) {
		exponent--;
	}
	int exponentExponent = (exponent==0) ? 1 : static_cast<int>(log10(abs(exponent))) + 1;
	if (abs(exponentExponent) > abs(exponentDigit)) {
		exponentDigit = abs(exponentExponent);
	}
	double mantissa = num / (pow(10, exponent));

	ss.setf(ios::showpos);
	ss.precision(precision+1);
	ss << fixed << mantissa << "e";
	ss.flags(MyUtil::IOMANIP_FLAGS_SAVED);
	ss << ( (exponent>=0) ? "+" : "-" );
	ss << setfill('0') << setw(exponentDigit) << right << abs(exponent);
	string result = ss.str();
	if (num>=0 && isPos==0) {
		result.replace(0, 1, " ");
	}
	return result;
}

string GetFNotation(double num, int precision, int isPos) {
	stringstream ss;
	ss.setf(ios::showpos);
	ss.precision(precision);
	ss << fixed << num;
	string result = ss.str();
	if (num>=0 && isPos==0) {
		result.replace(0, 1, " ");
	}
	return result;
}

ios::fmtflags SetIomanipScientificPos(int precision) {
	ResetIomanip();
	cout.setf(ios::showpos | ios::scientific);
	// cout.precision(precision);
	cout.precision(precision+1);
	return cout.flags();
}

void ResetIomanip() {
	cout.flags(MyUtil::IOMANIP_FLAGS_SAVED);
}

void _PrintErr(const string& s) {
	p("Error : ");
	pn(s);
	pn("Please press Enter key to continue...");
	getchar();
}


}
#endif
