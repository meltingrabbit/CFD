/*
言語:
C++
OS:
Microsoft Windows 10 Home (64 bit)
コンパイル:
C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall x86_amd64 & cl /EHsc /O2 /nologo $file_name
日付:
2017/01/12
*/
#ifndef MY_GEOMETRY_CLASS_H
#define MY_GEOMETRY_CLASS_H

namespace MyGeometry {

using namespace std;

template <typename T> class Point;
template <typename T> class Domain;

template <typename T>
class Point {
private:
	T x;
	T y;
	T z;
public:
	Point();
	Point(T x, T y, T z=0);
	inline T& X();
	inline T& Y();
	inline T& Z();
};

template <typename T>
class Domain {
private:
	Point<T> lb;		// left-bottom
	Point<T> rt;		// right-top
public:
	Domain();
	Domain(Point<T> lb, Point<T> rt);
	Domain(T x1, T y1, T x2, T y2);
	inline T& X1();
	inline T& X2();
	inline T& Y1();
	inline T& Y2();
};


// ########################################
// グローバル演算子オーバーロード

template <typename T>
inline ostream& operator<<(ostream& os, Point<T>& a) {
	return os << a.X() << "\t" << a.Y() << "\t" << a.Z();
}
template <typename T>
inline ostream& operator<<(ostream& os, Domain<T>& a) {
	return os << a.X1() << "\t" << a.Y1() << "\t" << a.X2() << "\t" << a.Y2();
}



// ########################################
// class実装

template <typename T>
Point<T>::Point():
	x(), y(), z()
{
}

template <typename T>
Point<T>::Point(T x, T y, T z):
	x(x), y(y), z(z)
{
}

template <typename T>
inline T& Point<T>::X() {
	return x;
}

template <typename T>
inline T& Point<T>::Y() {
	return y;
}

template <typename T>
inline T& Point<T>::Z() {
	return z;
}


template <typename T>
Domain<T>::Domain():
	lb(), rt()
{
}

template <typename T>
Domain<T>::Domain(Point<T> lb, Point<T> rt):
	lb(lb), rt(rt)
{
}

template <typename T>
Domain<T>::Domain(T x1, T y1, T x2, T y2):
	lb(x1,y1), rt(x2,y2)
{
}

template <typename T>
inline T& Domain<T>::X1() {
	return lb.X();
}

template <typename T>
inline T& Domain<T>::X2() {
	return rt.X();
}

template <typename T>
inline T& Domain<T>::Y1() {
	return lb.Y();
}

template <typename T>
inline T& Domain<T>::Y2() {
	return rt.Y();
}

}

#endif
