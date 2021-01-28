#include<bits/stdc++.h>
#include<iostream>
#include<iomanip>
#include<complex>
#include<fftw3.h>
using namespace std;

complex<double> mul(complex<double> a,complex<double> b)
{
	double x1 = a.real();
	double y1 = a.imag();
	double x2 = b.real();
	double y2 = b.imag();
	double realpart = x1*x2 - y1*y2;
	double imagpart = x1*y2 + y1*x2;
	complex<double> res(realpart,imagpart);
	return res;
}

int main()
{
	complex<double> a(1,1);
	complex<double> b(2,1);
	complex<double> res = mul(a,b);
	cout<<res.real()<<" "<<res.imag()<<endl;
	
	//we will now typecast this into a fftw_complex
	fftw_complex t;
	t[0] = res.real();
	t[1] = res.imag();
	cout<<t[0]<<" "<<t[1]<<endl;
	return 0;
}
