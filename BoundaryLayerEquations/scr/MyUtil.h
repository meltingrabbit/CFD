#ifndef MYUTIL_H
#define MYUTIL_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>		// ss
#include <iomanip>		// マニピュレータ
#include <string>
#include <time.h>

namespace MyUtil {

using namespace std;
const std::ios::fmtflags IOMANIP_FLAGS_SAVED = std::cout.flags();
const double PI = M_PI;

// 汎用Print
template <typename T> void p(T a);
template <typename T> void pn(T a);
void _ps();
template <typename T, typename... ARGV> void _ps(T first, ARGV... argv);
template <typename... ARGV> void ps(ARGV... argv);
template <typename... ARGV> void psn(ARGV... argv);
// 汎用数学関数
// アルゴリズム
// 数値計算
// 行列計算Eigen
// Util
string GetStringTime();

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

// ########################################
// 汎用数学関数

// ########################################
// アルゴリズム

// ########################################
// 数値計算

// ########################################
// Util
string GetStringTime() {
	stringstream ss;
	int y, m, d, hor, min, sec;
	time_t now = time(NULL);
	struct tm *pnow = localtime(&now);

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

}
#endif
