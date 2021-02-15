#include<bits/stdc++.h>
#include<iostream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<complex>
#include<fftw3.h>
#include<math.h>
using namespace std;

int main()
{
	int N1 = 8;
	int N2 = 8;
	double pi = 3.1428;
	
	/* rho(x,y) will be given to us and we need to calculate phi(x,y) */
	double L = 2*pi;
	double dx = L/N1;
	double dy = L/N2;
	
	double *rho_x_y;
	fftw_complex *rho_kx_ky;
	rho_x_y = (double*)malloc(sizeof(double)*N1*N2);
	rho_kx_ky = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);
	for(int n1=0;n1<N1;n1++)
	{
		for(int n2=0;n2<N2;n2++)
		{
			rho_x_y[n1*N2+n2] = sin(n1*dx)*cos(n2*dy);
		}
	}
	/*we will now find the fourier transform of this rho function.*/
	fftw_plan p;
	p = fftw_plan_dft_r2c_2d(N1,N2,rho_x_y,rho_kx_ky,FFTW_ESTIMATE);
	fftw_execute(p);
	
	/*we will find phi_kx_ky by dividing rho_kx_ky by (kx^2 + ky^2).*/
	fftw_complex *phi_kx_ky = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);

	for(int n1=0;n1<N1;n1++)
	{
		for(int n2=0;n2<N2;n2++)
		{
		   /* we will first calculate kx and ky. if both kx and ky equals 0 then then phi_kx_ky 
		    becomes infinite in this case since kx*kx + ky*ky = 0. we assume phi_kx_ky=0 when 
		    both kx and ky are zero.*/
		   double kx = n1*dx;
		   double ky = n2*dy;
           if(kx==0 && ky==0)
           {
			   phi_kx_ky[n1*N2+n2][0] = 0.0;
			   phi_kx_ky[n1*N2+n2][1] = 0.0; 
		   }
           else
           {
			  double temp = kx*kx + ky*ky;
			  phi_kx_ky[n1*N2+n2][0] = rho_kx_ky[n1*N2+n2][0]/temp;
			  phi_kx_ky[n1*N2+n2][1] = rho_kx_ky[n1*N2+n2][1]/temp;
		   }
		}
	}
	
/* we will now take inverse fourier transform of phi_kx_ky to get phi_x_y.Note that while taking
inverse fourier transform we should also do normalisation.*/
      double *phi_x_y;
      phi_x_y = (double*)malloc(sizeof(double)*N1*N2);
      fftw_plan q;
      q = fftw_plan_dft_c2r_2d(N1,N2,phi_kx_ky,phi_x_y,FFTW_ESTIMATE);
      fftw_execute(q);
      /*Note that we have to normalise after inverse fourier transform.*/
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  phi_x_y[n1*N2+n2] = phi_x_y[n1*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	  
	  /* we will print phi_x_y */
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  cout<<n1<<" "<<n2<<" "<<phi_x_y[n1*N1+n2]<<endl;
		  }
	  }
	  
	 fftw_destroy_plan(q);
	 fftw_destroy_plan(p);
	 fftw_cleanup();
	 return 0;
}
