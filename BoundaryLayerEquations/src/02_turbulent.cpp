#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>		// ss
#include <fstream>		// ifstream, ofstream
#include <algorithm>	// maxなど
#include <chrono>		// 時間計測
#include <vector>
#include "MyUtil.h"

using namespace std;

inline int outputFile(const string& filename, const ostringstream& oss);
int MakeGnuplotFile();
void PrintBar(int num);
void PrintDbar(int num);

// #############################################################
// 初期条件
// #############################################################
const double U0 = 20;					// 流速
// 以下，標準大気表高度0mより
const double RHO = 1.1255;				// 密度
const double NU = 1.4607 * pow(10,-5);	// 動粘性係数
const double MU = RHO * NU;				// 粘性係数
const double L = 0.5;					// 板の長さ
const double RE = RHO * U0 * L / MU;	// レイノルズ数
const double NH = 2;
const double H = NH * L * 0.37 / pow(RE, 0.2);	// 厚さ．deltaのNH倍
const int NX = 50000;					// X分割数
const int NY = 1000;						// Y分割数
const double dx = 1.0 * L / NX;
const double dy = 1.0 * H / NY;
const double A_PLUS = 26;
const double K_SMALL = 0.4;
const double K_LARGE = 0.0168;
const double C_CP = 1.6;
const double U_DIF =U0;
const double C_WK = 0.25;
const double C_KLEB = 0.3;
// #############################################################
// 除算回数削減のため
const double MUdRHO = MU / RHO;
const double dRHO = 1/RHO;
const double ddx = 1/dx;
const double ddy = 1/dy;
const double dU0 = 1/U0;
const double dMU = 1/MU;
const double dA_PLUS = A_PLUS;
const string stime = MyUtil::GetStringTime();
const string outputConditionFileName = "task_2_" + stime + "_Condition.txt";
const string outputThickness4CfFileName = "task_2_" + stime + "_Thickness4Cf.txt";
const string outputFlowVectorFileName = "task_2_" + stime + "_flowVector.txt";
const string outputFlowVectorAllFileName = "task_2_" + stime + "_flowVectorAll.txt";
const string outputErrFileName = "task_2_" + stime + "_Err.txt";
const string outputMutFileName = "task_2_" + stime + "_Mut.txt";
const string outputPlusFileName = "task_2_" + stime + "_Plus.txt";
const string outputGpFileName = "task_2_" + stime + ".gp";
// 一つ前も保存．2行行列
vector<vector<double>> u = vector<vector<double>>(2, vector<double>(NY+1, U0));
vector<vector<double>> v = vector<vector<double>>(2, vector<double>(NY+1, 0));
vector<double> mu_t_inn = vector<double>(NY+1, 0);
vector<double> mu_t_out = vector<double>(NY+1, 0);
vector<double> mu_t = vector<double>(NY+1, 0);
vector<double> dudy = vector<double>(NY+1, 0);				// j+1,jの差分．dudy[NY] = dudy[NY-1] とする

int main() {
	ios::sync_with_stdio(false);	// iostreamとstandard C streamsの同期をOFF

	// 初期化
	u[0][0] = 0;
	u[1][0] = 0;
	v[0][0] = 0;
	v[1][0] = 0;
	double bThickness, dThickness, mThickness, eThickness;
	std::chrono::system_clock::time_point startTime, endTime;
	double elapsedTime;

	cout.setf(ios::scientific);
	cout.precision(4);
	PrintDbar(79);
	MyUtil::pn("  TURBULENT BOUNDARY LAYER SIMULATION");
	PrintDbar(79);
	stringstream ss;
	MyUtil::psn("  ", "U0\t: ", U0, "\t[m/s]\t\t主流速度");
	MyUtil::psn("  ", "RE\t: ", RE, "\t[]\t\tレイノルズ数");
	MyUtil::psn("  ", "RHO\t: ",RHO,"\t[kg/m^3]\t密度");
	MyUtil::psn("  ", "NU\t: ", NU, "\t[m^2/s]\t\t粘性係数");
	MyUtil::psn("  ", "MU\t: ", MU, "\t[Pa s]\t\t動粘性係数");
	MyUtil::psn("  ", "L\t: ",  L,  "\t[m]\t\t幅");
	ss.str("");
	ss << NH;
	MyUtil::psn("  ", "H\t: ",  H,  "\t[m]", "\t\t厚さ (0.37 x L / RE^(0.2) x ", ss.str(), ")");
	ss.str("");
	ss << setw(10) << right << NX;
	MyUtil::psn("  ", "NX\t: ", ss.str(), "\t[]\t\t幅分割数");
	ss.str("");
	ss << setw(10) << right << NY;
	MyUtil::psn("  ", "NY\t: ", ss.str(), "\t[]\t\t厚分割数");
	MyUtil::psn("  ", "dx\t: ", dx, "\t[m]\t\t格子幅");
	MyUtil::psn("  ", "dy\t: ", dy, "\t[m]\t\t格子厚");
	MyUtil::psn("  ", "stime\t: ", stime, "\t\tプログラム起動時間（ログID）");
	PrintDbar(79);
	std::cout.flags(MyUtil::IOMANIP_FLAGS_SAVED);
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：計算開始");
	startTime = std::chrono::system_clock::now();

	// ファイル出力バッファー
	ostringstream osFlowVector, osFlowVectorAll, osThickness4Cf, osMut, osPlus, osErr;

	cout.setf(ios::showpos | ios::scientific);
	cout.precision(4);
	bThickness = 0;
	for (int i=0; i<=NX; ++i) {		// i=0はgiven
		// 積分の初期化
		dThickness = 0;
		mThickness = 0;
		eThickness = 0;

		// #############################
		// u, v の計算
		// #############################
		u[1][0] = 0;
		v[1][0] = 0;
		if (i != 0) {
			double sumV = 0;
			// j : 0 ～ NY-1
			for (int j=1; j<NY; ++j) {		// j=0,NYはgiven
				// 配列アクセスを減らす
				// Previous, Current, Next
				double upp = u[0][j-1];
				double upc = u[0][j];
				double upn = u[0][j+1];
				double vpc = v[0][j];
				double ucp = u[1][j-1];
				double ucc;
				double mup = MU + mu_t[j-1];
				double muc = MU + mu_t[j  ];
				double mun = MU + mu_t[j+1];
				ucc = upc + (dx/upc) * ( - (vpc * (upn - upp) * ddy / 2) );
				// 粘性項
				ucc += (dx/upc) * dRHO * ddy * (ddy/2) * ( (muc+mun)*(upn-upc) - (mup+muc)*(upc-upp) );
				sumV += ( -((ucp-upp)*ddx) - ((ucc-upc)*ddx) );
				u[1][j] = ucc;
				v[1][j] = sumV * dy / 2;
				// 境界層厚さの探索
				if (u[1][j-1] < U0*0.995 && ucc >= U0*0.995) {
					// 線形補間
					int idxL = j - 1;
					int idxR = j;
					double idxDbl = idxL + (U0*0.995-u[1][idxL])/(u[1][idxR]-u[1][idxL]) * (idxR - idxL);
					bThickness = idxDbl * dy;
				}
				// 排除，運動量，エネルギー厚さの計算
				// dyは後でかける
				double udU0 = ucc * dU0;
				dThickness += (1-udU0);
				mThickness += udU0*(1-udU0);
				eThickness += udU0*(1-udU0*udU0);
				// dudyの計算
				dudy[j-1] = (u[1][j]-u[1][j-1]) * ddy;
			}
			// j : 0 ～ NY
			double upp = u[0][NY-1];
			double upc = u[0][NY];
			double ucp = u[1][NY-1];
			double ucc = U0;
			sumV += ( -((ucp-upp)*ddx) - ((ucc-upc)*ddx) );
			u[1][NY] = ucc;
			v[1][NY] = sumV * dy / 2;
			// dudyの計算
			dudy[NY-1] = (u[1][NY]-u[1][NY-1]) * ddy;
			dudy[NY] = dudy[NY-1];			// 境界なので計算できない
		}

		if (i%5000 == 0) {
			for (int j=0; j<=NY; ++j) {
				osFlowVectorAll << i*dx << "\t" << j*dy << "\t" << u[1][j] << "\t" << v[1][j] << "\n";
				if (j%50 == 0) {
					osFlowVector << i*dx << "\t" << j*dy << "\t" << u[1][j] << "\t" << v[1][j] << "\n";
				}
			}
			osFlowVector << "\n\n";
			osFlowVectorAll << "\n\n";
		}

		// 表面摩擦係数の計算
		double tauW = MU * u[1][1] * ddy;
		double cf = tauW * 2 * dRHO * dU0*dU0;

		osThickness4Cf << i*dx << "\t" << bThickness << "\t" << dThickness*dy << "\t" << mThickness*dy << "\t" << eThickness*dy << "\t" << cf << "\n";

		// 理論値の境界層厚さとの誤差の計算
		double bThicknessTheory = 0.37*i*dx / pow(RHO * U0 * i*dx * dMU, 0.2);
		osErr << i*dx << "\t" << abs( (bThickness - bThicknessTheory)/bThicknessTheory );
		// 理論値の壁面摩擦係数との誤差の計算
		double cfTheory = 0.0592 / pow(RHO * U0 * i*dx * dMU, 0.2);
		osErr << "\t" << abs( (cf - cfTheory)/cfTheory ) << "\n";

		// #############################
		// μ の計算
		// #############################
		double yPlusdy = sqrt(RHO * dMU * dudy[0]);
		// fMax, yMaxの探査
		double fMax = 0;
		int yMaxIdx = 0;				// yMax index
		for (int j=0; j<=NY; ++j) {
			double y = j*dy;
			double f = y * abs(dudy[j]) * (1 - exp(-(yPlusdy * y * dA_PLUS)));
			if (fMax < f) {
				fMax = f;
				yMaxIdx = j;
			}
		}
		double yMax = yMaxIdx * dy;

		int isMutInnOverMutOut = 0;
		for (int j=0; j<=NY; ++j) {
			// y_crossを見つけ次第，mu_t_innは計算しない
			double y = j*dy;
			// mu_t_out の計算
			double fWake = min(yMax*fMax, C_WK*yMax*U_DIF*U_DIF/fMax);
			double fKleb = 1.0 / (1 + 5.5 * pow(C_KLEB*y/yMax,6));
			mu_t_out[j] = K_LARGE * C_CP * RHO * fWake * fKleb;
			if (isMutInnOverMutOut == 0) {
				// mu_t_inn の計算
				double l = K_SMALL * y * (1 - exp(-(yPlusdy * y * dA_PLUS)));
				mu_t_inn[j] = RHO * l*l * abs(dudy[j]);
				if (mu_t_out[j] < mu_t_inn[j]) {
					mu_t[j] = mu_t_out[j];
					isMutInnOverMutOut = 1;
				} else {
					mu_t[j] = mu_t_inn[j];
				}
			} else {
				mu_t[j] = mu_t_out[j];
			}
		}

		if (i%1000 == 0) {
			for (int j=0; j<=NY; ++j) {
				osMut << i*dx << "\t" << j*dy << "\t" << mu_t_inn[j] << "\t" << mu_t_out[j] << "\n";
			}
			osMut << "\n\n";
		}

		// 現在の行を一つ前へシフト
		u[0].swap(u[1]);
		v[0].swap(v[1]);
	}

	// x=Lでの u+, y+, mu_tのファイル出力
	// u,vはすでにswap済みなので，[0]を使う
	double uTau = sqrt((MU * dudy[0]) / RHO);
	for (int j=0; j<=NY; ++j) {
		double uPlus = u[0][j] / uTau;
		double y = j*dy;
		double yPlus = sqrt(RHO * dMU * dudy[0]) * y;
		osPlus << y << "\t" << yPlus << "\t" << uPlus << "\t" << mu_t[j] * dMU << "\n";
	}

	std::cout.flags(MyUtil::IOMANIP_FLAGS_SAVED);

	endTime = std::chrono::system_clock::now();
	elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count();
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：計算終了    処理時間 ", elapsedTime, " [ms]");
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：出力開始");
	startTime = std::chrono::system_clock::now();

	// ファイル出力
	outputFile(outputThickness4CfFileName, osThickness4Cf);
	outputFile(outputFlowVectorFileName, osFlowVector);
	outputFile(outputFlowVectorAllFileName, osFlowVectorAll);
	outputFile(outputPlusFileName, osPlus);
	outputFile(outputMutFileName, osMut);
	outputFile(outputErrFileName, osErr);
	ostringstream osCondition;
	osCondition << "U0\t" << U0 << "\n";
	osCondition << "RHO\t" << RHO << "\n";
	osCondition << "NU\t" << NU << "\n";
	osCondition << "MU\t" << MU << "\n";
	osCondition << "L\t" << L << "\n";
	osCondition << "RE\t" << RE << "\n";
	osCondition << "NH\t" << NH << "\n";
	osCondition << "H\t" << H << "\n";
	osCondition << "NX\t" << NX << "\n";
	osCondition << "NY\t" << NY << "\n";
	osCondition << "dx\t" << dx << "\n";
	osCondition << "dy\t" << dy << "\n";
	outputFile(outputConditionFileName, osCondition);
	MakeGnuplotFile();

	endTime = std::chrono::system_clock::now();
	elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count();
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：出力終了    処理時間 ", elapsedTime, " [ms]");
	return 0;
}

inline int outputFile(const string& filename, const ostringstream& oss) {
	fstream fs;
	fs.open(filename, ios::out);
	if(! fs.is_open()) {
		return EXIT_FAILURE;
	}
	fs << oss.str();
	fs.close();
	return 1;
}

int MakeGnuplotFile() {
	// ###########################################################
	// 可視化のためのgnuplot load fileを作成する関数．
	// プログラムの主要素ではないので省略．
	// ###########################################################
}

void PrintDbar(int num) {
	string dbar(num, '=');
	MyUtil::pn(dbar);
}
