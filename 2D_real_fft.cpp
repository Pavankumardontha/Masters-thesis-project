#include<bits/stdc++.h>
#include<complex>
#include<iomanip>
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
using namespace std;
/* function is sin(3x)*sin(5y). We are multiplying here */
int main()
{
    int N1 = 8;
    int N2 = 8;
    double pi = 3.1428;
    /*f(x,y) is always real.Its fourier transform may be complex.*/
    double L = 2*pi;
    /*we will now find dx and dy.*/
    double dx = L/N1;
    double dy = L/N2;
    
    double *u_x_y;//this defines our function
    fftw_complex *u_kx_ky;
    u_x_y = (double*)malloc(sizeof(double)*N1*N2);
    u_kx_ky = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);
    for(int n1=0;n1<N1;n1++)
    {
		for(int n2=0;n2<N2;n2++)
		{
			u_x_y[n1*N2+n2] = sin(3*dx*n1) + sin(2*dy*n2);
		}
	}
	/* we will define a plan and execute it.*/
	fftw_plan p;
	p = fftw_plan_dft_r2c_2d(N1,N2,u_x_y,u_kx_ky,FFTW_ESTIMATE);
	fftw_execute(p);
	
	/*lets print the result.*/
	for(int n1=0;n1<N1;n1++)
	 {
		 for(int n2=0;n2<N2;n2++)
		 {
/*value of n1 and n2 followed by real part and imaginary part of u_kx_ky followed by magnitude of u_kx_ky.*/ 
			 cout<<setw(12)<<n1<<" "<<n2<<"\t"<<setw(12)<<u_kx_ky[n1*N2+n2][0]<<"\t"<<setw(12)<<u_kx_ky[n1*N2+n2][1]<<"\t"<<setw(12)<<sqrt(pow(u_kx_ky[n1*N2+n2][0],2)+pow(u_kx_ky[n1*N2+n2][1],2))<<endl;
		 }
	 }
	 
	 /* now we will try to do inverse fourier transform and get back our original function u_x_y
	  from u_kx_ky.*/
	  double *ui_x_y;
	  ui_x_y = (double*)malloc(sizeof(double)*N1*N2);
	  /* we will define inverse plan and execute it.*/
	  fftw_plan q;
	  q = fftw_plan_dft_c2r_2d(N1,N2,u_kx_ky,ui_x_y,FFTW_ESTIMATE);
	  fftw_execute(q);
	  
	  /* we will print the result.*/
	  cout<<"Un-Normalised Inverse Fourier Transform of u_x_y"<<endl;
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  cout<<n1<<" "<<n2<<" "<<u_x_y[n1*N2+n2]<<" "<<ui_x_y[n1*N2+n2]<<endl;
		  }
	  }
	  
	  /*Note that we have to normalise after inverse fourier transform.*/
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  ui_x_y[n1*N2+n2] = ui_x_y[n1*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	  cout<<"Normalised Inverse Fourier Transform of u_x_y"<<endl;
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  cout<<n1<<" "<<n2<<" "<<u_x_y[n1*N2+n2]<<" "<<ui_x_y[n1*N2+n2]<<endl;
		  }
	  }
	 /*We can clearly see that after normalisation the values of ui_x_y and u_x_y are same.*/
	  
	 fftw_destroy_plan(q);
	 fftw_destroy_plan(p);
	 fftw_cleanup();
	 return 0;
	
}
