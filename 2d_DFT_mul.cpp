#include<bits/stdc++.h>
#include<complex>
#include<iomanip>
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
using namespace std;
/* function is sin(3x)*sin(5x). We are multiplying here */
int main()
{
     int N1 = 8;
     int N2 = 8;
     /*
     function f(x,y)=sin(((2*pie)/N1)*3*n1) * sin(((2*pie)/N2)*5*n2); 
     n1 = [0,..N1-1],n2=[0,1,..N2-1].
     */
     fftw_complex *in;
     fftw_complex *out;
     in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);
     out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);
	 fftw_plan p;
	 for(int n1=0;n1<N1;n1++)
	 {
		 for(int n2=0;n2<N2;n2++)
		 {
			 in[n1*N2+n2][0] = sin(((2*M_PI)/N1)*3*n1) * sin(((2*M_PI)/N2)*5*n2);
			 in[n1*N2+n2][1] = 0.0;
		 }
	 }	
	 p = fftw_plan_dft_2d(N1,N2,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	 fftw_execute(p);
	 for(int n1=0;n1<N1;n1++)
	 {
		 for(int n2=0;n2<N2;n2++)
		 {
			 cout<<n1<<" "<<n2<<"\t"<<setw(12)<<out[n1*N2+n2][0]<<"\t"<<setw(12)<<out[n1*N2+n2][1]<<"\t"<<setw(12)<<sqrt(pow(out[n1*N2+n2][0],2)+pow(out[n1*N2+n2][1],2))<<endl;
		 }
	 }
	 fftw_destroy_plan(p);
	 fftw_cleanup();
	 return 0;
}
