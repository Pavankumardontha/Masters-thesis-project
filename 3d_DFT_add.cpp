#include<bits/stdc++.h>
#include<iostream>
#include<stdlib.h>
#include<fftw3.h>
#include<stdio.h>
#include<complex>
#include<iomanip>
using namespace std;

int main()
{
	int N1=5;
	int N2=5;
	int N3=5;
	/*
     function f(x,y)=sin(((2*pie)/N1)*3*n1) * sin(((2*pie)/N2)*5*n2) * sin(((2*pie)/N3)*2*n3); 
     n1 = [0,..N1-1],n2=[0,1,..N2-1],n3=[0,1...N3-1].
    */
    fftw_complex *in;
    fftw_complex *out;
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
    fftw_plan p;
    for(int n1=0;n1<N1;n1++)
    {
		for(int n2=0;n2<N2;n2++)
		{
			for(int n3=0;n3<N3;n3++)
			{
				in[n1*n2*N3+n3][0] = sin(((2*M_PI)/N1)*2*n1) + sin(((2*M_PI)/N2)*3*n2) + sin(((2*M_PI)/N3)*4*n3);
			    in[n1*n2*N3+n3][1] = 0.0;	
			}
		}
	}
	p = fftw_plan_dft_3d(N1,N2,N3,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for(int n1=0;n1<N1;n1++)
	{
		for(int n2=0;n2<N2;n2++)
		{
			for(int n3=0;n3<N3;n3++)
			{
				cout<<n1<<" "<<n2<<" "<<n3<<"\t"<<setw(12)<<out[n1*n2*N3+n3][0]<<"\t"<<setw(12)<<out[n1*n2*N3+n3][1]<<"\t"<<setw(12)<<sqrt(pow(out[n1*n2*N3+n3][0],2)+pow(out[n1*n2*N3+n3][1],2))<<endl;
			}
		}
	}
	fftw_destroy_plan(p);
	fftw_cleanup();
	return 0;
	
}
