/*
言語:
C++
OS:
Microsoft Windows 10 Home (64 bit)
コンパイル:
C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall x86_amd64 & cl /EHsc /O2 /nologo $file_name
日付:
2016/11/21
*/
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <sstream>		// ss
#include <algorithm>	// maxなど
#include "MyUtil.h"
// #include "MyRange.h"
#include "MyGeometry.h"
// #include "MyTime.h"

// #define NO_PRINT

using namespace std;
using namespace MyUtil;
using namespace MyGeometry;

void SetIndex();
void SetPBoundaryCondition();
void SetUVBoundaryCondition();
void SolvePoissonEquation(int tStep);
void SolveVelocityEquation(int tStep);
void CalcCoefficient(int tStep, double time);
double XIdx2X(int idx);
double YIdx2Y(int idx);


const double RE    = 70.0;			// Re
const double invRE = 1.0/RE;
const double CFL   = 0.2;			// CFL Number クーラン数
Domain<double> domCalcArea(-10.0,-10.0,30.0,10.0);
Domain<double> domBox(-0.5,-0.5,0.5,0.5);
const double U0   = 1.0;			// 流速
const double V0   = 0.0;			// 流速
const double P0   = 0.0;			// 圧力
const int T_STEP_END = 50000;
const int PRINT_CONSOLE_MOD = 5;
const int OUTPUT_FILE_MOD = 50;

const double dx = 0.1;
const double dy = 0.1;
double dt;
const double invDx = 1.0/dx;
const double invDy = 1.0/dy;
double invDt;
// index
Domain<int> IDX_CALC_AREA;
Domain<int> IDX_BOX;
Point<int>  IDX_CENTER;
const int   IDX_MARGIN = 2;			// 差分法のためのマージン
int NX;
int NY;
vector<vector<double>> U;
vector<vector<double>> V;
vector<vector<double>> P;

// 出力用にやむを得ずグローバル
int ITR_P;
double RES_P;

const string outputCoefficientFileName = "task_1_Coefficient.dat";


int main() {
	ios::sync_with_stdio(false);	// iostreamとstandard C streamsの同期をOFF

	SetIndex();			// 境界のindex計算
	// pn(IDX_CENTER);
	// pn(IDX_CALC_AREA);
	// pn(IDX_BOX);
	// pn(XIdx2X(100));
	// pn(XIdx2X(200));
	// pn(NX);
	// pn(NY);


	// クーラン条件
	dt = CFL * min(dx,dy);
	invDt = 1.0/dt;

	// 初期化
	U = vector<vector<double>>(NX+1, vector<double>(NY+1, U0));
	V = vector<vector<double>>(NX+1, vector<double>(NY+1, V0));
	P = vector<vector<double>>(NX+1, vector<double>(NY+1, P0));
	outputFile(outputCoefficientFileName, "");

	pn("t STEP\tt\tITR_P\tRES_P\t\tC_D\tC_L\tC_P_1\tC_P_2");
	pn(GetDbar(72));


	// // ノイズ入れとく？
	// {
	// 	int i1 = IDX_BOX.X1();
	// 	int i2 = IDX_BOX.X2();
	// 	int j1 = IDX_BOX.Y1();
	// 	int j2 = IDX_BOX.Y2();
	// 	REP(k,30) {
	// 		P[i1+k][j2+2] = 0.5;
	// 	}
	// }

	// return 1;

	double time = 0;
	SetPBoundaryCondition();
	SetUVBoundaryCondition();

	string filenameU = "task_1_2dVec_U_tStep=" + ToString(0) + ".csv";
	output2dVector2Csv(filenameU, U);
	string filenameV = "task_1_2dVec_V_tStep=" + ToString(0) + ".csv";
	output2dVector2Csv(filenameV, V);
	string filenameP = "task_1_2dVec_P_tStep=" + ToString(0) + ".csv";
	output2dVector2Csv(filenameP, P);

	// REP(tStep,T_STEP_END+1) {
	for (int tStep=1;tStep<=T_STEP_END;++tStep) {
		time += dt;
		#ifndef NO_PRINT
			if (tStep%PRINT_CONSOLE_MOD == 0) {
				ps(tStep, "\t", GetFNotation(time,1,1), "\t");
			}
		#endif
		SolvePoissonEquation(tStep);
		// SetPBoundaryCondition();
		SolveVelocityEquation(tStep);
		// SetUVBoundaryCondition();
		CalcCoefficient(tStep, time);
		#ifndef NO_PRINT
			if (tStep%PRINT_CONSOLE_MOD == 0) {
				pn("");
			}
		#endif

		if (tStep%OUTPUT_FILE_MOD == 0) {
			string filenameU = "task_1_2dVec_U_tStep=" + ToString(tStep) + ".csv";
			output2dVector2Csv(filenameU, U);
			string filenameV = "task_1_2dVec_V_tStep=" + ToString(tStep) + ".csv";
			output2dVector2Csv(filenameV, V);
			string filenameP = "task_1_2dVec_P_tStep=" + ToString(tStep) + ".csv";
			output2dVector2Csv(filenameP, P);
			// string filename;
			// filename = "task_1_P_tStep=" + ToString(tStep) + ".dat";
			// outputFile(filename, Matrix2GnuplotPm3dMapData(P, dx, dy));
		}
	}


// 	// MyTime sw;
// 	cin >> dt;
// 	sw.End();
// 	pn(sw.Watch<MyTime::ns>());
// 	pn(sw.Watch<MyTime::us>());
// 	pn(sw.Watch<MyTime::ms>());
// 	pn(sw.Watch<MyTime::s>());
// 	pn(sw.Watch<MyTime::m>());
// 	pn(sw.Watch<MyTime::h>());
// // #include <typeinfo>
// 	auto sss = sw.Watch<MyTime::ns>();
// 	// pn(typeid(sss).name());

	return 0;
}


void SetIndex() {
	int nx,ny;

	IDX_CALC_AREA.X1() = 0;
	IDX_CALC_AREA.Y1() = 0;

	nx = static_cast<int>(-domCalcArea.X1()/dx);
	ny = static_cast<int>(-domCalcArea.Y1()/dy);
	IDX_CENTER.X() = IDX_CALC_AREA.X1() + nx;
	IDX_CENTER.Y() = IDX_CALC_AREA.Y1() + ny;

	nx = static_cast<int>(domCalcArea.X2()/dx);
	ny = static_cast<int>(domCalcArea.Y2()/dy);
	IDX_CALC_AREA.X2() = IDX_CENTER.X() + nx;
	IDX_CALC_AREA.Y2() = IDX_CENTER.Y() + ny;

	nx = static_cast<int>(-domBox.X1()/dx);
	ny = static_cast<int>(-domBox.Y1()/dy);
	IDX_BOX.X1() = IDX_CENTER.X() - nx;
	IDX_BOX.Y1() = IDX_CENTER.Y() - ny;

	nx = static_cast<int>(domBox.X2()/dx);
	ny = static_cast<int>(domBox.Y2()/dy);
	IDX_BOX.X2() = IDX_CENTER.X() + nx;
	IDX_BOX.Y2() = IDX_CENTER.Y() + ny;

	NX = IDX_CALC_AREA.X2();
	NY = IDX_CALC_AREA.Y2();
}

void SetPBoundaryCondition() {
	static int i1 = IDX_BOX.X1();
	static int i2 = IDX_BOX.X2();
	static int j1 = IDX_BOX.Y1();
	static int j2 = IDX_BOX.Y2();

	// 境界
	REP(i,NX+1) {
		P[i][ 0] = P0;
		P[i][ 1] = P0;
		P[i][NY  ] = P0;
		P[i][NY-1] = P0;
	}
	REP(j,NY+1) {
		P[ 0][j] = P0;
		P[ 1][j] = P0;
		P[NX  ][j] = P0;
		P[NX-1][j] = P0;
	}

	// psn(i1,"\t",i2,"\t",j1,"\t",j2,"\t",NX,"\t",NY,"\t");
// return;
	// 物体境界
	P[i1][j1] = P[i1-1][j1-1];
	P[i1][j2] = P[i1-1][j2+1];
	P[i2][j1] = P[i2+1][j1-1];
	P[i2][j2] = P[i2+1][j2+1];

	for (int j=j1+1;j<j2;++j) {
		P[i1][j] = P[i1-1][j];
		P[i2][j] = P[i2+1][j];
	}
	for (int i=i1+1;i<i2;++i) {
		P[i][j1] = P[i][j1-1];
		P[i][j2] = P[i][j2+1];
	}
}

void SetUVBoundaryCondition() {
	static int i1 = IDX_BOX.X1();
	static int i2 = IDX_BOX.X2();
	static int j1 = IDX_BOX.Y1();
	static int j2 = IDX_BOX.Y2();

	// 境界
	REP(j,NY+1) {
		U[0][j] = U0;
		V[0][j] = V0;
		U[1][j] = U0;
		V[1][j] = V0;
		U[NX-1][j] = 2.0 * U[NX-2][j] - U[NX-3][j];
		V[NX-1][j] = 2.0 * V[NX-2][j] - V[NX-3][j];
		U[NX  ][j] = 2.0 * U[NX-1][j] - U[NX-2][j];
		V[NX  ][j] = 2.0 * V[NX-1][j] - V[NX-2][j];
	}
	REP(i,NX+1) {
		U[i][ 1]   = 2.0 * U[i][   2] - U[i][   3];
		V[i][ 1]   = 2.0 * V[i][   2] - V[i][   3];
		U[i][ 0]   = 2.0 * U[i][   1] - U[i][   2];
		V[i][ 0]   = 2.0 * V[i][   1] - V[i][   2];
		U[i][NY-1] = 2.0 * U[i][NY-2] - U[i][NY-3];
		V[i][NY-1] = 2.0 * V[i][NY-2] - V[i][NY-3];
		U[i][NY  ] = 2.0 * U[i][NY-1] - U[i][NY-2];
		V[i][NY  ] = 2.0 * V[i][NY-1] - V[i][NY-2];
	}
// return;
	// 物体境界
	for (int j=j1;j<=j2;++j) {
		U[i1][j] = 0.0;
		U[i2][j] = 0.0;
		V[i1][j] = 0.0;
		V[i2][j] = 0.0;
	}
	for (int i=i1;i<=i2;++i) {
		U[i][j1] = 0.0;
		U[i][j2] = 0.0;
		V[i][j1] = 0.0;
		V[i][j2] = 0.0;
	}
	// U[i1+1][j1+1] = -U[i1-1][j1-1];
	// U[i1+1][j2-1] = -U[i1-1][j2+1];
	// U[i2-1][j1+1] = -U[i2+1][j1-1];
	// U[i2-1][j2-1] = -U[i2+1][j2+1];
	// V[i1+1][j1+1] = -V[i1-1][j1-1];
	// V[i1+1][j2-1] = -V[i1-1][j2+1];
	// V[i2-1][j1+1] = -V[i2+1][j1-1];
	// V[i2-1][j2-1] = -V[i2+1][j2+1];
	// for (int j=j1+2;j<=j2-2;++j) {
	// 	U[i1+1][j] = -U[i1-1][j];
	// 	U[i2-1][j] = -U[i2+1][j];
	// 	V[i1+1][j] = -V[i1-1][j];
	// 	V[i2-1][j] = -V[i2+1][j];
	// }
	// for (int i=i1+2;i<=i2-2;++i) {
	// 	U[i][j1+1] = -U[i][j1-1];
	// 	U[i][j2-1] = -U[i][j2+1];
	// 	V[i][j1+1] = -V[i][j1-1];
	// 	V[i][j2-1] = -V[i][j2+1];
	// }
}

void SolvePoissonEquation(int tStep) {
	// 外側1格子は解けない
	// static vector<vector<double>> rhs = vector<vector<double>>(NX+1, vector<double>(NY+1, 0));
	vector<vector<double>> rhs = vector<vector<double>>(NX+1, vector<double>(NY+1, 0));
	static int MAX_ITERATION = 100;
	static int i1 = IDX_BOX.X1();
	static int i2 = IDX_BOX.X2();
	static int j1 = IDX_BOX.Y1();
	static int j2 = IDX_BOX.Y2();
	double ux,uy,vx,vy;

	// ループ判定がめんどいので，不要な物体内部も計算 ?? どうしよ．計算コスト？
	for (int i=0+2;i<=NX-2;++i) {
		for (int j=0+2;j<=NY-2;++j) {
			if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) continue;
			// if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) {psn(i,"\t",j);continue;}
			// if ((i > i1 && i < i2) && (j > j1 && j < j2)) continue;
			ux = (U[i+1][j] - U[i-1][j]) *0.5 * invDx;
			uy = (U[i][j+1] - U[i][j-1]) *0.5 * invDy;
			vx = (V[i+1][j] - V[i-1][j]) *0.5 * invDx;
			vy = (V[i][j+1] - V[i][j-1]) *0.5 * invDy;
			rhs[i][j] = (ux + vy) *invDt - (ux*ux + 2.0 * uy * vx + vy*vy);
			// rhs[i][j] = (ux + vy) *invDt - (ux * ux + 2.0 * ux * uy + uy * uy);
			// ux = (U[i+1][j] - U[i-1][j]) / (2.0 * dx);
			// uy = (U[i][j+1] - U[i][j-1]) / (2.0 * dy);
			// vx = (V[i+1][j] - V[i-1][j]) / (2.0 * dx);
			// vy = (V[i][j+1] - V[i][j-1]) / (2.0 * dy);
			// rhs[i][j] = (ux + vy) / dt - (ux * ux + 2.0 * ux * uy + uy * uy);
		}
	}

	double dp;
	double res = 0.0;		// 残差
	int itrNum = MAX_ITERATION+1;		// breakしなかったときの回数
	static const double OMEGA = 1.00;			// ポアソン方程式緩和係数
	static const double EPSILON = 0.0001;
	// invese (2/dx/dx + 2/dy/dy)
	static const double invDenom = 1.0 / (2.0*invDx*invDx + 2.0*invDy*invDy);
	REP(itr,MAX_ITERATION) {
		res = 0.0;
		for (int i=0+2;i<=NX-2;++i) {
			for (int j=0+2;j<=NY-2;++j) {
				if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) continue;
				// if ((i > i1 && i < i2) && (j > j1 && j < j2)) continue;

				dp = (P[i+1][j] + P[i-1][j]) * invDx*invDx + (P[i][j+1] + P[i][j-1]) * invDy*invDy - rhs[i][j];
				// dp = (P[i+1][j] + P[i-1][j]) / pow(dx,2.0) + (P[i][j+1] + P[i][j-1]) / pow(dy,2.0) - rhs[i][j];
				dp = dp * invDenom - P[i][j];
				// dp = dp / (2.0 / pow(dx, 2.0) + 2.0 / pow(dy, 2.0)) - P[i][j];
				res += dp*dp;
				// pn(dp);
				P[i][j] += OMEGA * dp;
			}
		}

		// 圧力境界のアップデート
		SetPBoundaryCondition();
		// psn((2.0*invDx*invDx * 2.0*invDy*invDy), "\t", 2/(dx*dx)+2/(dy*dy));
		// ps(res,"\t");
		res = sqrt(res / static_cast<double>(NX*NY));
		// ps(res, "\t");
		// psn(NX*NY);
		if (res < EPSILON) {
			itrNum = itr+1;
			break;
		}
	}

	ITR_P = itrNum;
	RES_P = res;
	#ifndef NO_PRINT
		if (tStep%PRINT_CONSOLE_MOD == 0) {
			ps(itrNum, "\t", GetENotation(res,2,1,2), "\t");
		}
	#endif
	return;
}

void SolveVelocityEquation(int tStep) {
	// 外側2格子は解けない
	// static vector<vector<double>> rhsU = vector<vector<double>>(NX+1, vector<double>(NY+1, 0));
	// static vector<vector<double>> rhsV = vector<vector<double>>(NX+1, vector<double>(NY+1, 0));
	vector<vector<double>> rhsU = vector<vector<double>>(NX+1, vector<double>(NY+1, 0));
	vector<vector<double>> rhsV = vector<vector<double>>(NX+1, vector<double>(NY+1, 0));
	static int i1 = IDX_BOX.X1();
	static int i2 = IDX_BOX.X2();
	static int j1 = IDX_BOX.Y1();
	static int j2 = IDX_BOX.Y2();
	static double inv12 = 1.0 / 12.0;
	static double inv4 = 1.0 / 4.0;

	for (int i=0+2;i<=NX-2;++i) {
		for (int j=0+2;j<=NY-2;++j) {
			if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) continue;
			rhsU[i][j] = -(P[i+1][j] - P[i-1][j]) * 0.5 * invDx;
			rhsV[i][j] = -(P[i][j+1] - P[i][j-1]) * 0.5 * invDy;
		}
	}

	for (int i=0+2;i<=NX-2;++i) {
		for (int j=0+2;j<=NY-2;++j) {
			if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) continue;
			rhsU[i][j] += (U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) * invRE * invDx*invDx;
			rhsU[i][j] += (U[i][j+1] - 2.0 * U[i][j] + U[i][j-1]) * invRE * invDy*invDy;
			rhsV[i][j] += (V[i+1][j] - 2.0 * V[i][j] + V[i-1][j]) * invRE * invDx*invDx;
			rhsV[i][j] += (V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) * invRE * invDy*invDy;
		}
	}

	// KawamuraScheme用，内部線形外挿x方向）
	for (int j=j1+1;j<=j2-1;++j) {
		U[i1+1][j] = 2.0 * U[i1][j] - U[i1-1][j];
		U[i2-1][j] = 2.0 * U[i2][j] - U[i2+1][j];
		V[i1+1][j] = 2.0 * V[i1][j] - V[i1-1][j];
		V[i2-1][j] = 2.0 * V[i2][j] - V[i2+1][j];
	}

	for (int i=0+2;i<=NX-2;++i) {
		for (int j=0+2;j<=NY-2;++j) {
			if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) continue;
			rhsU[i][j] +=  -U[i][j] * (-U[i+2][j] + 8.0 * (U[i+1][j] - U[i-1][j]) + U[i-2][j]) * inv12 * invDx;
			rhsU[i][j] += -abs(U[i][j]) * (U[i+2][j] - 4.0 * U[i+1][j] + 6.0 * U[i][j] - 4.0 * U[i-1][j] + U[i-2][j]) * inv4 * invDx;
			rhsV[i][j] +=  -U[i][j] * (-V[i+2][j] + 8.0 * (V[i+1][j] - V[i-1][j]) + V[i-2][j]) * inv12 * invDx;
			rhsV[i][j] += -abs(U[i][j]) * (V[i+2][j] - 4.0 * V[i+1][j] + 6.0 * V[i][j] - 4.0 * V[i-1][j] + V[i-2][j]) * inv4 * invDx;
		}
	}

	// KawamuraScheme用，内部線形外挿y方向）
	for (int i=i1+1;i<=i2-1;++i) {
		U[i][j1+1] = 2.0 * U[i][j1] - U[i][j1-1];
		U[i][j2-1] = 2.0 * U[i][j2] - U[i][j2+1];
		V[i][j1+1] = 2.0 * V[i][j1] - V[i][j1-1];
		V[i][j2-1] = 2.0 * V[i][j2] - V[i][j2+1];
	}

	for (int i=0+2;i<=NX-2;++i) {
		for (int j=0+2;j<=NY-2;++j) {
			if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) continue;
			rhsU[i][j] +=  -V[i][j] * (-U[i][j+2] + 8.0 * (U[i][j+1] - U[i][j-1]) + U[i][j-2]) * inv12 * invDy;
			rhsU[i][j] += -abs(V[i][j]) * (U[i][j+2] - 4.0 * U[i][j+1] + 6.0 * U[i][j] - 4.0 * U[i][j-1] + U[i][j-2]) * inv4 * invDy;
			rhsV[i][j] +=  -V[i][j] * (-V[i][j+2] + 8.0 * (V[i][j+1] - V[i][j-1]) + V[i][j-2]) * inv12 * invDy;
			rhsV[i][j] += -abs(V[i][j]) * (V[i][j+2] - 4.0 * V[i][j+1] + 6.0 * V[i][j] - 4.0 * V[i][j-1] + V[i][j-2]) * inv4 * invDy;
		}
	}

	for (int i=0+2;i<=NX-2;++i) {
		for (int j=0+2;j<=NY-2;++j) {
			// 物体内部と外側1列
			if ((i1 <= i && i <= i2) && (j1 <= j && j <= j2)) continue;
			U[i][j] += dt * rhsU[i][j];
			V[i][j] += dt * rhsV[i][j];
		}
	}


	// 速度境界のアップデート
	SetUVBoundaryCondition();
}

void CalcCoefficient(int tStep, double time) {
	static int i1 = IDX_BOX.X1();
	static int i2 = IDX_BOX.X2();
	static int j1 = IDX_BOX.Y1();
	static int j2 = IDX_BOX.Y2();

	double cp1,cp2,cd,cl;
	ostringstream ossCoeff;

	// Cp = 2*p

	cd = 0.0;
	for (int j=j1;j<j2;++j) {
		double cpFront = P[i1][j] + P[i1][j+1];
		double cpBack  = P[i2][j] + P[i2][j+1];
		cd += (cpFront - cpBack) * dy;
	}

	cl = 0.0;
	for (int i=i1+1;i<i2;++i) {
		double cpBtm = P[i][j1] + P[i+1][j1];
		double cpTop = P[i][j2] + P[i+1][j2];
		cl += (cpBtm - cpTop) * dx;
	}

	cp1 = 2.0 * P[i2+i2-i1][j1];
	cp2 = 2.0 * P[i2+i2-i1][j2];

	#ifndef NO_PRINT
		if (tStep%PRINT_CONSOLE_MOD == 0) {
			ps(GetFNotation(cd,3,1),"\t",GetFNotation(cl,3,1),"\t",GetFNotation(cp1,3,1),"\t",GetFNotation(cp2,3,1));
			// ps(P[i1][j1],"\t",P[i2][j2]);
		}
	#endif

	ossCoeff << tStep << "\t" << time << "\t" << ITR_P << "\t" << GetENotation(RES_P,2,1,2) << "\t" << cd << "\t" << cl << "\t" << cp1 << "\t" << cp2 << "\n";
	outputFile(outputCoefficientFileName, ossCoeff, 1);
}


double XIdx2X(int idx) {
	return dx*(idx-IDX_CENTER.X());
}

double YIdx2Y(int idx) {
	return dy*(idx-IDX_CENTER.Y());
}

