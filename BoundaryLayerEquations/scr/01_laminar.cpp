#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>		// ss
#include <fstream>		// ifstream, ofstream
#include <chrono>		// 時間計測
#include <vector>
#include "MyUtil.h"

using namespace std;

inline double GetVWall(double x);
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
const double H = NH * L * 5.3 / sqrt(RE);	// 厚さ．ブラジウス厚さのNH倍
const int NX = 50000;					// X分割数
// const int NY = 250;						// Y分割数 1.1節
const int NY = 200;						// Y分割数 1.2節
// const int NY = 200;						// Y分割数 1.3節
const double dx = 1.0 * L / NX;
const double dy = 1.0 * H / NY;
const double V_WALL = 0;						// 境界層吸い込み 1.1節
// const double V_WALL = 0.012;						// 境界層吸い込み 1.2節
// const double V_WALL = -0.012;						// 境界層吸い込み 1.2節
// const double V_WALL = 0;						// 境界層吸い込み 1.3節
const double V_WALL_RANGE_L = 0.0;				// 1.1節
const double V_WALL_RANGE_R = 0.5;				// 1.1節
// const double V_WALL_RANGE_L = 0.2;				// 1.2節
// const double V_WALL_RANGE_R = 0.4;				// 1.2節
// const double V_WALL_RANGE_L = 0.0;				// 1.3節
// const double V_WALL_RANGE_R = 0.5;				// 1.3節
const double dpdx = 0;						// 圧力勾配 1.1節
// const double dpdx = 0;						// 圧力勾配 1.2節
// const double dpdx = 95;						// 圧力勾配 1.3節
// #############################################################
// 除算回数削減のため
const double MUdRHO = MU / RHO;
const double dRHO = 1/RHO;
const double dMU = 1/MU;
const double ddx = 1/dx;
const double ddy = 1/dy;
const double dU0 = 1/U0;
const string stime = MyUtil::GetStringTime();
const string outputConditionFileName = "task_1_" + stime + "_Condition.txt";
const string outputErrFileName = "task_1_" + stime + "_Err.txt";
const string outputThickness4CfFileName = "task_1_" + stime + "_Thickness4Cf.txt";
const string outputFlowVectorFileName = "task_1_" + stime + "_flowVector.txt";
const string outputFlowVectorAllFileName = "task_1_" + stime + "_flowVectorAll.txt";
const string outputGpFileName = "task_1_" + stime + ".gp";
// 一つ前も保存．2行行列
vector<vector<double>> u = vector<vector<double>>(2, vector<double>(NY+1, U0));
vector<vector<double>> v = vector<vector<double>>(2, vector<double>(NY+1, 0));

int main() {
	ios::sync_with_stdio(false);	// iostreamとstandard C streamsの同期をOFF

	// 初期化
	u[0][0] = 0;
	u[1][0] = 0;
	v[0][0] = GetVWall(0);
	v[1][0] = GetVWall(0);
	double bThickness, dThickness, mThickness, eThickness;
	std::chrono::system_clock::time_point startTime, endTime;
	double elapsedTime;

	cout.setf(ios::scientific);
	cout.precision(4);
	PrintDbar(79);
	MyUtil::pn("  LAMINAR BOUNDARY LAYER SIMULATION");
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
	MyUtil::psn("  ", "H\t: ",  H,  "\t[m]", "\t\t厚さ (Blasius BLT x ", ss.str(), ")");
	ss.str("");
	ss << setw(10) << right << NX;
	MyUtil::psn("  ", "NX\t: ", ss.str(), "\t[]\t\t幅分割数");
	ss.str("");
	ss << setw(10) << right << NY;
	MyUtil::psn("  ", "NY\t: ", ss.str(), "\t[]\t\t厚分割数");
	MyUtil::psn("  ", "dx\t: ", dx, "\t[m]\t\t格子幅");
	MyUtil::psn("  ", "dy\t: ", dy, "\t[m]\t\t格子厚");
	MyUtil::ps("  ", "V_WALL: ", V_WALL, "\t[m/s]");
	std::cout.flags(MyUtil::IOMANIP_FLAGS_SAVED);
	MyUtil::psn("\t\t境界層吸い込み (", V_WALL_RANGE_L, " <= x <= ", V_WALL_RANGE_R, ")");
	cout.setf(ios::scientific);
	cout.precision(4);
	MyUtil::psn("  ", "dpdx\t: ", dpdx, "\t[Pa/m]\t\t圧力勾配");
	MyUtil::psn("  ", "stime\t: ", stime, "\t\tプログラム起動時間（ログID）");
	PrintDbar(79);
	std::cout.flags(MyUtil::IOMANIP_FLAGS_SAVED);
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：計算開始");
	startTime = std::chrono::system_clock::now();

	// ファイル出力バッファー
	ostringstream osFlowVector, osFlowVectorAll, osThickness4Cf, osErr;

	cout.setf(ios::showpos | ios::scientific);
	cout.precision(4);
	bThickness = 0;
	for (int i=0; i<=NX; ++i) {		// i=0はgiven
		// 積分の初期化
		dThickness = 0;
		mThickness = 0;
		eThickness = 0;

		u[1][0] = 0;
		double vWall = GetVWall(i*dx);
		v[1][0] = vWall;
		if (i == 0) {
			u[0][0] = U0;
			u[1][0] = U0;
		}
		if (i != 0) {
			double sumV = 0;
			double mainstreamU = (dpdx == 0 ? U0 : u[0][NY]);			// 圧力勾配の有無で場合分け
			double dmainstreamU = 1/mainstreamU;
			for (int j=1; j<NY; ++j) {		// j=0,NYはgiven
				// 配列アクセスを減らす
				// Previous, Current, Next
				double upp = u[0][j-1];
				double upc = u[0][j];
				double upn = u[0][j+1];
				double vpc = v[0][j];
				double ucp = u[1][j-1];
				double ucc;
				ucc = upc + (dx/upc) * ( (MUdRHO * (upn - 2*upc + upp) * ddy*ddy) - (vpc * (upn - upp) * ddy / 2) );
				// 圧力勾配を付加
				ucc -= (dx/upc) * ( dRHO * dpdx );
				sumV += ( -((ucp-upp)*ddx) - ((ucc-upc)*ddx) );
				u[1][j] = ucc;
				v[1][j] = vWall + sumV * dy / 2;
				// 境界層厚さの探索
				if (ucp < mainstreamU*0.995 && ucc >= mainstreamU*0.995) {
					// 線形補間
					int idxL = j - 1;
					int idxR = j;
					double idxDbl = idxL + (mainstreamU*0.995-u[1][idxL])/(u[1][idxR]-u[1][idxL]) * (idxR - idxL);
					bThickness = idxDbl * dy;
				}
				// 排除，運動量，エネルギー厚さの計算
				// dyは後でかける
				double udU0 = ucc * dmainstreamU;
				dThickness += (1-udU0);
				mThickness += udU0*(1-udU0);
				eThickness += udU0*(1-udU0*udU0);
			}
			double upp = u[0][NY-1];
			double upc = u[0][NY];
			double ucp = u[1][NY-1];
			double ucc = (dpdx == 0 ? U0 : ucp);			// 圧力勾配の有無で境界条件を変更
			sumV += ( -((ucp-upp)*ddx) - ((ucc-upc)*ddx) );
			u[1][NY] = ucc;
			v[1][NY] = vWall + sumV * dy / 2;
		}

		if (i%5000 == 0) {
			for (int j=0; j<=NY; ++j) {
				osFlowVectorAll << i*dx << "\t" << j*dy << "\t" << u[1][j] << "\t" << v[1][j] << "\n";
				if (j%10 == 0) {
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

		// Blasiusの境界層厚さとの誤差の計算
		double bThicknessTheory = 5.3*sqrt(i*dx) / sqrt(RHO * U0 * dMU);
		osErr << i*dx << "\t" << abs( (bThickness - bThicknessTheory)/bThicknessTheory );
		// Blasiusの壁面摩擦係数との誤差の計算
		double cfTheory = 0.664 / sqrt(RHO * U0 * i*dx * dMU);
		osErr << "\t" << abs( (cf - cfTheory)/cfTheory ) << "\n";

		// 現在の行を一つ前へシフト
		u[0].swap(u[1]);
		v[0].swap(v[1]);
	}
	std::cout.flags(MyUtil::IOMANIP_FLAGS_SAVED);

	endTime = std::chrono::system_clock::now();
	elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count();
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：計算終了    処理時間 ", elapsedTime, " [ms]");
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：出力開始");
	startTime = std::chrono::system_clock::now();

	// ファイル出力
	outputFile(outputErrFileName, osErr);
	outputFile(outputThickness4CfFileName, osThickness4Cf);
	outputFile(outputFlowVectorFileName, osFlowVector);
	outputFile(outputFlowVectorAllFileName, osFlowVectorAll);
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
	osCondition << "V_WALL\t" << V_WALL << "\n";
	osCondition << "V_WALL_RANGE_L\t" << V_WALL_RANGE_L << "\n";
	osCondition << "V_WALL_RANGE_R\t" << V_WALL_RANGE_R << "\n";
	osCondition << "dpdx\t" << dpdx << "\n";
	outputFile(outputConditionFileName, osCondition);
	MakeGnuplotFile();

	endTime = std::chrono::system_clock::now();
	elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count();
	MyUtil::psn("  ", MyUtil::GetStringTime(), "：出力終了    処理時間 ", elapsedTime, " [ms]");
	return 0;
}

inline double GetVWall(double x) {
	return (V_WALL_RANGE_L <= x && x <= V_WALL_RANGE_R) ? V_WALL : 0;
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
