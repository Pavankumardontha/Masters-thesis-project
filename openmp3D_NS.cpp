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
	int N1 = 50;
	int N2 = 50;
	int N3 = 50;
	double pi = 3.1428;
	double L = 2*pi;
	double dx = L/N1;
	double dy = L/N2;
	double dz = L/N3;
	double kx;
	double ky;
	double kz;
	double temp;
	int Nt = 20;
	double time;
	double time_max=1.0;
	double time_min=0;
	double dt=0.1;
	double v = 0.01; //some random value for viscosity
	
	int val = fftw_init_threads();
    fftw_plan_with_nthreads(13);
	
	/* every velocity component is a function of (x,y,z).let us define them */
	
	//present components
	double *ux = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *uy = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *uz = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	//old components
	double *ux_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *uy_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *uz_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	double *dx_ux_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dy_ux_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dz_ux_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	double *dx_uy_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dy_uy_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dz_uy_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	double *dx_uz_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dy_uz_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dz_uz_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	double *dxx_ux_o =(double*)malloc(sizeof(double)*N1*N2*N3);
	double *dyy_ux_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dzz_ux_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	double *dxx_uy_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dyy_uy_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dzz_uy_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	double *dxx_uz_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dyy_uz_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *dzz_uz_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	//we dont need fourier transforms of x1,x2 and x3 */
	double *x1_o = (double*)malloc(sizeof(double)*N1*N2*N3); 
	double *x2_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	double *x3_o = (double*)malloc(sizeof(double)*N1*N2*N3);
	
	/* We will consider initial conditions of taylor green vortex given by 
	 ux_o(x,y,z) = (cosx)*(siny)*(cosz)
	 uy_o(x,y,z) = -(sinx)*(cosy)*(cosz)
	 uz_o(x,y,z) = 0
	* /
	
	/*initialise the functions u_x_o , u_y_o and u_z_o.*/
	#pragma omp parallel for private(n1) private(n2)
	for(int n1=0;n1<N1;n1++)
	{
		for(int n2=0;n2<N2;n2++)
		{
			for(int n3=0;n3<N3;n3++)
			{
				ux_o[n1*n2*N3 + n3] = cos(n1*dx)*sin(n2*dy)*cos(n3*dz);
				ux_o[n1*n2*N3 + n3] = sin(n1*dx)*cos(n2*dy)*cos(n3*dz);
				uz_o[n1*n2*N3 + n3] = 0.0;
			}
		}
	}
	
	/* we will now define and allocate space to the fourier transforms of u_x , u_y and u_z */
	fftw_complex *f_ux_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3); //fourier transform of u_x_o
	fftw_complex *f_uy_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3); //fourier transform of u_y_o
	fftw_complex *f_uz_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3); //fourier transform of u_z_o
	
	fftw_complex *f_dx_ux_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3); 
	fftw_complex *f_dy_ux_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dz_ux_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	
	fftw_complex *f_dx_uy_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3); 
	fftw_complex *f_dy_uy_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dz_uy_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	
	fftw_complex *f_dx_uz_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3); 
	fftw_complex *f_dy_uz_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dz_uz_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	
	fftw_complex *f_dxx_ux_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dyy_ux_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dzz_ux_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
		
	fftw_complex *f_dxx_uy_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dyy_uy_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dzz_uy_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	
	fftw_complex *f_dxx_uz_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dyy_uz_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	fftw_complex *f_dzz_uz_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2*N3);
	
	#pragma omp parallel for
	for(time=time_min+dt;time<=time_max;time=time+dt)
	{
		/*we have to calculate total 18 terms.Lets calculate the terms in order so that we dont 
		 take confuse or miss any term. We will calculate first all the necessary terms required 
		 to calculate x1 then x2 and x3 respectively.So we require 18 plans. */
		  fftw_plan a1;
	      a1 = fftw_plan_dft_r2c_3d(N1,N2,N3,ux_o,f_ux_o,FFTW_ESTIMATE);
	      fftw_execute(a1);
	      fftw_plan a2;
	      a2 = fftw_plan_dft_r2c_3d(N1,N2,N3,uy_o,f_uy_o,FFTW_ESTIMATE);
	      fftw_execute(a2);
	      fftw_plan a3;
	      a3 = fftw_plan_dft_r2c_3d(N1,N2,N3,uz_o,f_uz_o,FFTW_ESTIMATE);
	      fftw_execute(a3);
	      /* we have all the fourier transforms of ux_o,uy_o and uz_o. We will now try to find 
	       other terms. */
	       
	      #pragma omp parallel for private(n1) private(n2)
	      for(int n1=0;n1<N1;n1++)
	      {
			  for(int n2=0;n2<N2;n2++)
			  {
				  for(int n3=0;n3<N3;n3++)
				  {
					  kx = n1*dx;
					  ky = n2*dy;
					  kz = n3*dz;
					  
					  f_dxx_ux_o[n1*n2*N3+n3][0] = -(kx*kx)*(f_ux_o[n1*n2*N3+n3][0]);
					  f_dxx_ux_o[n1*n2*N3+n3][1] = -(kx*kx)*(f_ux_o[n1*n2*N3+n3][1]);
					  f_dyy_ux_o[n1*n2*N3+n3][0] = -(ky*ky)*(f_ux_o[n1*n2*N3+n3][0]);
					  f_dyy_ux_o[n1*n2*N3+n3][1] = -(ky*ky)*(f_ux_o[n1*n2*N3+n3][1]);
					  f_dzz_ux_o[n1*n2*N3+n3][0] = -(kz*kz)*(f_ux_o[n1*n2*N3+n3][0]);
					  f_dzz_ux_o[n1*n2*N3+n3][1] = -(kz*kz)*(f_ux_o[n1*n2*N3+n3][1]);
					  f_dxx_uy_o[n1*n2*N3+n3][0] = -(kx*kx)*(f_uy_o[n1*n2*N3+n3][0]);
					  f_dxx_uy_o[n1*n2*N3+n3][1] = -(kx*kx)*(f_uy_o[n1*n2*N3+n3][1]);
					  f_dyy_uy_o[n1*n2*N3+n3][0] = -(ky*ky)*(f_uy_o[n1*n2*N3+n3][0]);
					  f_dyy_uy_o[n1*n2*N3+n3][1] = -(ky*ky)*(f_uy_o[n1*n2*N3+n3][1]);
					  f_dzz_uy_o[n1*n2*N3+n3][0] = -(kz*kz)*(f_uy_o[n1*n2*N3+n3][0]);
					  f_dzz_uy_o[n1*n2*N3+n3][1] = -(kz*kz)*(f_uy_o[n1*n2*N3+n3][1]);
					  f_dxx_uz_o[n1*n2*N3+n3][0] = -(kx*kx)*(f_uz_o[n1*n2*N3+n3][0]);
					  f_dxx_uz_o[n1*n2*N3+n3][1] = -(kx*kx)*(f_uz_o[n1*n2*N3+n3][1]);
					  f_dyy_uz_o[n1*n2*N3+n3][0] = -(ky*ky)*(f_uz_o[n1*n2*N3+n3][0]);
					  f_dyy_uz_o[n1*n2*N3+n3][1] = -(ky*ky)*(f_uz_o[n1*n2*N3+n3][1]);
					  f_dzz_uz_o[n1*n2*N3+n3][0] = -(kz*kz)*(f_uz_o[n1*n2*N3+n3][0]);
					  f_dzz_uz_o[n1*n2*N3+n3][1] = -(kz*kz)*(f_uz_o[n1*n2*N3+n3][1]);
					  /* now will calculate the single derivative terms*/
					  
					  /* derivatives of ux*/
					  f_dx_ux_o[n1*n2*N3+n3][0] = kx*(f_ux_o[n1*n2*N3+n3][0]);
					  f_dx_ux_o[n1*n2*N3+n3][1] = kx*(f_ux_o[n1*n2*N3+n3][1]);
					  temp = f_dx_ux_o[n1*n2*N3+n3][0];
					  f_dx_ux_o[n1*n2*N3+n3][0] = -f_dx_ux_o[n1*n2*N3+n3][1];
					  f_dx_ux_o[n1*n2*N3+n3][1] = temp;
					  
					  f_dy_ux_o[n1*n2*N3+n3][0] = ky*(f_ux_o[n1*n2*N3+n3][0]);
					  f_dy_ux_o[n1*n2*N3+n3][1] = ky*(f_ux_o[n1*n2*N3+n3][1]);
					  temp = f_dy_ux_o[n1*n2*N3+n3][0];
					  f_dy_ux_o[n1*n2*N3+n3][0] = -f_dy_ux_o[n1*n2*N3+n3][1];
					  f_dy_ux_o[n1*n2*N3+n3][1] = temp;
					  
					  f_dz_ux_o[n1*n2*N3+n3][0] = kz*(f_ux_o[n1*n2*N3+n3][0]);
					  f_dz_ux_o[n1*n2*N3+n3][1] = kz*(f_ux_o[n1*n2*N3+n3][1]);
					  temp = f_dz_ux_o[n1*n2*N3+n3][0];
					  f_dz_ux_o[n1*n2*N3+n3][0] = -f_dz_ux_o[n1*n2*N3+n3][1];
					  f_dz_ux_o[n1*n2*N3+n3][1] = temp;
					  
					  /* derivatives of uy */
					  f_dx_uy_o[n1*n2*N3+n3][0] = kx*(f_uy_o[n1*n2*N3+n3][0]);
					  f_dx_uy_o[n1*n2*N3+n3][1] = kx*(f_uy_o[n1*n2*N3+n3][1]);
					  temp = f_dx_uy_o[n1*n2*N3+n3][0];
					  f_dx_uy_o[n1*n2*N3+n3][0] = -f_dx_uy_o[n1*n2*N3+n3][1];
					  f_dx_uy_o[n1*n2*N3+n3][1] = temp;
					  
					  f_dy_uy_o[n1*n2*N3+n3][0] = ky*(f_uy_o[n1*n2*N3+n3][0]);
					  f_dy_uy_o[n1*n2*N3+n3][1] = ky*(f_uy_o[n1*n2*N3+n3][1]);
					  temp = f_dy_uy_o[n1*n2*N3+n3][0];
					  f_dy_uy_o[n1*n2*N3+n3][0] = -f_dy_uy_o[n1*n2*N3+n3][1];
					  f_dy_uy_o[n1*n2*N3+n3][1] = temp;
					  
					  f_dz_uy_o[n1*n2*N3+n3][0] = kz*(f_uy_o[n1*n2*N3+n3][0]);
					  f_dz_uy_o[n1*n2*N3+n3][1] = kz*(f_uy_o[n1*n2*N3+n3][1]);
					  temp = f_dz_uy_o[n1*n2*N3+n3][0];
					  f_dz_uy_o[n1*n2*N3+n3][0] = -f_dz_uy_o[n1*n2*N3+n3][1];
					  f_dz_uy_o[n1*n2*N3+n3][1] = temp;
					  
					  /* derivatives of uz*/
					  f_dx_uz_o[n1*n2*N3+n3][0] = kx*(f_uz_o[n1*n2*N3+n3][0]);
					  f_dx_uz_o[n1*n2*N3+n3][1] = kx*(f_uz_o[n1*n2*N3+n3][1]);
					  temp = f_dx_uz_o[n1*n2*N3+n3][0];
					  f_dx_uz_o[n1*n2*N3+n3][0] = -f_dx_uz_o[n1*n2*N3+n3][1];
					  f_dx_uz_o[n1*n2*N3+n3][1] = temp;
					  
					  f_dy_uz_o[n1*n2*N3+n3][0] = ky*(f_uz_o[n1*n2*N3+n3][0]);
					  f_dy_uz_o[n1*n2*N3+n3][1] = ky*(f_uz_o[n1*n2*N3+n3][1]);
					  temp = f_dy_uz_o[n1*n2*N3+n3][0];
					  f_dy_uz_o[n1*n2*N3+n3][0] = -f_dy_uz_o[n1*n2*N3+n3][1];
					  f_dy_uz_o[n1*n2*N3+n3][1] = temp;
					  
					  f_dz_uz_o[n1*n2*N3+n3][0] = kz*(f_uz_o[n1*n2*N3+n3][0]);
					  f_dz_uz_o[n1*n2*N3+n3][1] = kz*(f_uz_o[n1*n2*N3+n3][1]);
					  temp = f_dz_uz_o[n1*n2*N3+n3][0];
					  f_dz_uz_o[n1*n2*N3+n3][0] = -f_dz_uz_o[n1*n2*N3+n3][1];
					  f_dz_uz_o[n1*n2*N3+n3][1] = temp;
				  }
			  }
		  }
		  
	      /* we have to now take the inverse fourier transforms */
	      fftw_plan a4 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dxx_ux_o,dxx_ux_o,FFTW_ESTIMATE);
	      fftw_plan a5 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dyy_ux_o,dyy_ux_o,FFTW_ESTIMATE);
	      fftw_plan a6 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dzz_ux_o,dzz_ux_o,FFTW_ESTIMATE);
	      
	      fftw_plan a7 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dxx_uy_o,dxx_ux_o,FFTW_ESTIMATE);
	      fftw_plan a8 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dyy_uy_o,dyy_ux_o,FFTW_ESTIMATE);
	      fftw_plan a9 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dzz_uy_o,dzz_ux_o,FFTW_ESTIMATE);
	      
	      fftw_plan a10 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dxx_uz_o,dxx_ux_o,FFTW_ESTIMATE);
	      fftw_plan a11 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dyy_uz_o,dyy_ux_o,FFTW_ESTIMATE);
	      fftw_plan a12 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dzz_uz_o,dzz_ux_o,FFTW_ESTIMATE);
	      
	      fftw_plan a13 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dx_ux_o,dx_ux_o,FFTW_ESTIMATE);
	      fftw_plan a14 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dy_ux_o,dy_ux_o,FFTW_ESTIMATE);
	      fftw_plan a15 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dz_ux_o,dz_ux_o,FFTW_ESTIMATE);
	      
	      fftw_plan a16 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dx_uy_o,dx_uy_o,FFTW_ESTIMATE);
	      fftw_plan a17 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dy_uy_o,dy_uy_o,FFTW_ESTIMATE);
	      fftw_plan a18 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dz_uy_o,dz_uy_o,FFTW_ESTIMATE);
	      
	      fftw_plan a19 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dx_uz_o,dx_uz_o,FFTW_ESTIMATE);
	      fftw_plan a20 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dy_uz_o,dy_uz_o,FFTW_ESTIMATE);
	      fftw_plan a21 = fftw_plan_dft_c2r_3d(N1,N2,N3,f_dz_uz_o,dz_uz_o,FFTW_ESTIMATE);
		
		  /* we have to now execute all these plans one by one. */
		  fftw_execute(a4);
		  fftw_execute(a5);
		  fftw_execute(a6);
		  
		  fftw_execute(a7);
		  fftw_execute(a8);
		  fftw_execute(a9);
		  
		  fftw_execute(a10);
		  fftw_execute(a11);
		  fftw_execute(a12);
		  
		  fftw_execute(a13);
		  fftw_execute(a14);
		  fftw_execute(a15);
		  
		  fftw_execute(a16);
		  fftw_execute(a17);
		  fftw_execute(a18);
		  
		  fftw_execute(a19);
		  fftw_execute(a20);
		  fftw_execute(a21);
		  
	    /*normalising all the terms obtained by inverse fourier transforms. */
	    #pragma omp parallel for private(n1) private(n2)
		for(int n1=0;n1<N1;n1++)
	    {
		  for(int n2=0;n2<N2;n2++)
		  {
			  for(int n3=0;n3<N3;n3++)
			  {
				  dx_ux_o[n1*n2*N3 + n3] = dx_ux_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dy_ux_o[n1*n2*N3 + n3] = dy_ux_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dz_ux_o[n1*n2*N3 + n3] = dz_ux_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dxx_ux_o[n1*n2*N3 + n3] = dxx_ux_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dyy_ux_o[n1*n2*N3 + n3] = dyy_ux_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dzz_ux_o[n1*n2*N3 + n3] = dzz_ux_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  
				  dx_uy_o[n1*n2*N3 + n3] = dx_uy_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dy_uy_o[n1*n2*N3 + n3] = dy_uy_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dz_uy_o[n1*n2*N3 + n3] = dz_uy_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dxx_uy_o[n1*n2*N3 + n3] = dxx_uy_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dyy_uy_o[n1*n2*N3 + n3] = dyy_uy_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dzz_uy_o[n1*n2*N3 + n3] = dzz_uy_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  
				  dx_uz_o[n1*n2*N3 + n3] = dx_uz_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dy_uz_o[n1*n2*N3 + n3] = dy_uz_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dz_uz_o[n1*n2*N3 + n3] = dz_uz_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dxx_uz_o[n1*n2*N3 + n3] = dxx_uz_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dyy_uz_o[n1*n2*N3 + n3] = dyy_uz_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
				  dzz_uz_o[n1*n2*N3 + n3] = dzz_uz_o[n1*n2*N3 + n3]/(float(N1)*float(N2)*float(N3));
			  }
		  }
	  }
		
		/* calculating x1_o , x2_o and x3_o */
		#pragma omp parallel for private(n1) private(n2)
		for(int n1=0;n1<N1;n1++)
		{
			for(int n2=0;n2<N2;n2++)
			{
				for(int n3=0;n3<N3;n3++)
				{
					x1_o[n1*n2*N3+n3] = v*(dxx_ux_o[n1*n2*N3+n3]+dyy_ux_o[n1*n2*N3+n3]+dzz_ux_o[n1*n2*N3+n3]);
					x1_o[n1*n2*N3+n3] = x1_o[n1*n2*N3+n3] - (ux_o[n1*n2*N3+n3]*dx_ux_o[n1*n2*N3+n3]);
					x1_o[n1*n2*N3+n3] = x1_o[n1*n2*N3+n3] - (uy_o[n1*n2*N3+n3]*dy_ux_o[n1*n2*N3+n3]);
					x1_o[n1*n2*N3+n3] = x1_o[n1*n2*N3+n3] - (uz_o[n1*n2*N3+n3]*dz_ux_o[n1*n2*N3+n3]);
					
					x2_o[n1*n2*N3+n3] = v*(dxx_uy_o[n1*n2*N3+n3]+dyy_uy_o[n1*n2*N3+n3]+dzz_uy_o[n1*n2*N3+n3]);
					x2_o[n1*n2*N3+n3] = x2_o[n1*n2*N3+n3] - (ux_o[n1*n2*N3+n3]*dx_uy_o[n1*n2*N3+n3]);
					x2_o[n1*n2*N3+n3] = x2_o[n1*n2*N3+n3] - (uy_o[n1*n2*N3+n3]*dy_uy_o[n1*n2*N3+n3]);
					x2_o[n1*n2*N3+n3] = x2_o[n1*n2*N3+n3] - (uz_o[n1*n2*N3+n3]*dz_uy_o[n1*n2*N3+n3]);
					
					x3_o[n1*n2*N3+n3] = v*(dxx_uz_o[n1*n2*N3+n3]+dyy_uz_o[n1*n2*N3+n3]+dzz_uz_o[n1*n2*N3+n3]);
					x3_o[n1*n2*N3+n3] = x3_o[n1*n2*N3+n3] - (ux_o[n1*n2*N3+n3]*dx_uz_o[n1*n2*N3+n3]);
					x3_o[n1*n2*N3+n3] = x3_o[n1*n2*N3+n3] - (uy_o[n1*n2*N3+n3]*dy_uz_o[n1*n2*N3+n3]);
					x3_o[n1*n2*N3+n3] = x3_o[n1*n2*N3+n3] - (uz_o[n1*n2*N3+n3]*dz_uz_o[n1*n2*N3+n3]);
					
				}
			}
		}
		
		/*calculating the velocitys at present time */
		#pragma omp parallel for private(n1) private(n2)
		for(int n1=0;n1<N1;n1++)
		{
			for(int n2=0;n2<N2;n2++)
			{
				for(int n3=0;n3<N3;n3++)
				{
					ux[n1*n2*N3+n3] = ux_o[n1*n2*N3+n3] + ((x1_o[n1*n2*N3+n3])*dt);
					uy[n1*n2*N3+n3] = uy_o[n1*n2*N3+n3] + ((x2_o[n1*n2*N3+n3])*dt);
					uz[n1*n2*N3+n3] = uz_o[n1*n2*N3+n3] + ((x3_o[n1*n2*N3+n3])*dt);
					//we can copy all the components into some file here 
				}
			}
		}
		
		/* we will print the results.*/
		cout<<"time:"<<time<<endl;
		for(int n1=0;n1<N1;n1++)
		{
			for(int n2=0;n2<N2;n2++)
			{
				for(int n3=0;n3<N3;n3++)
				{
					cout<<n1<<" "<<n2<<" "<<n3<<" "<<ux[n1*n2*N3+n3]<<endl;
				}
			}
		}
		
		
		
		/*assign the present components to the old one's */ 
		#pragma omp parallel for private(n1) private(n2)
		for(int n1=0;n1<N1;n1++)
		{
			for(int n2=0;n2<N2;n2++)
			{
				for(int n3=0;n3<N3;n3++)
				{
					ux_o[n1*n2*N3+n3] = ux[n1*n2*N3+n3];
					uy_o[n1*n2*N3+n3] = uy[n1*n2*N3+n3];
					uz_o[n1*n2*N3+n3] = uz[n1*n2*N3+n3];
				}
			}
		}
		
		/* we have to destroy all the plans which we have created till now.*/
		fftw_destroy_plan(a1);
		fftw_destroy_plan(a2);
		fftw_destroy_plan(a3);
		
		fftw_destroy_plan(a4);
		fftw_destroy_plan(a5);
		fftw_destroy_plan(a6);
		
		fftw_destroy_plan(a7);
		fftw_destroy_plan(a8);
		fftw_destroy_plan(a9);
		
		fftw_destroy_plan(a10);
		fftw_destroy_plan(a11);
		fftw_destroy_plan(a12);
		
		fftw_destroy_plan(a13);
		fftw_destroy_plan(a14);
		fftw_destroy_plan(a15);
		
		fftw_destroy_plan(a16);
		fftw_destroy_plan(a17);
		fftw_destroy_plan(a18);
		
		fftw_destroy_plan(a19);
		fftw_destroy_plan(a20);
		fftw_destroy_plan(a21);
		}
		return 0;
}
