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
#ifndef MY_RANGE_CLASS_H
#define MY_RANGE_CLASS_H

// #include "MyUtil.h"
#include <utility>			// pair

using namespace std;

template <typename T>
class MyRange {
public:
	using range_t = pair<T, T>;
private:
	range_t range;
public:
	MyRange(T l, T r);
	inline T& L();			// left, begin
	inline T& R();			// right, end
};

template <typename T>
MyRange<T>::MyRange(T l, T r):
	range(l, r)
{
}

template <typename T>
inline T& MyRange<T>::L() {
	return range.first;
}

template <typename T>
inline T& MyRange<T>::R() {
	return range.second;
}


#endif
