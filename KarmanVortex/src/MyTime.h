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
#ifndef MY_TIME_CLASS_H
#define MY_TIME_CLASS_H

#include <chrono>		// 時間計測

using namespace std;

class MyTime {
public:
	using ns = chrono::nanoseconds;
	using us = chrono::microseconds;
	using ms = chrono::milliseconds;
	using s  = chrono::seconds;
	using m  = chrono::minutes;
	using h  = chrono::hours;
private:
	chrono::system_clock::time_point startTime;
	chrono::system_clock::time_point endTime;
public:
	MyTime();
	void Start();
	void End();
	template <typename T> auto Elapsed();
	template <typename T> auto Watch();
};

MyTime::MyTime()
{
	startTime = chrono::system_clock::now();
}

void MyTime::Start() {
	startTime = chrono::system_clock::now();
}

void  MyTime::End() {
	endTime = chrono::system_clock::now();
}

template <typename T>
auto MyTime::Elapsed() {
	return chrono::duration_cast<T>(endTime-startTime).count();
}

template <typename T>
auto MyTime::Watch() {
	End();
	return Elapsed<T>();
}

#endif
