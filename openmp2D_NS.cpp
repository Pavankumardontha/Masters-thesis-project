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
	int N1 = 500;
	int N2 = 500;
	double pi = 3.1428;
	
	double L = 2*pi;
	double dx = L/N1;
	double dy = L/N2;
	int val = fftw_init_threads();
    fftw_plan_with_nthreads(3);
    
/*lets say omega is omega_z=sinx*cosy.Note that in 2D omege has only z-component and hence acts 
like a scalar.*/

	double* omega_x_y_o;
	double* omega_x_y;
	fftw_complex *omega_kx_ky_o;
	omega_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
	omega_x_y = (double*)malloc(sizeof(double)*N1*N2);
	omega_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);
	
	#pragma omp parallel for private(n2)
	for(int n1=0;n1<N1;n1++)
	{
		for(int n2=0;n2<N2;n2++)
		{
			omega_x_y_o[n1*N2+n2] = sin(n1*dx)*cos(n2*dy);
		}
	}
/*we need dx_omega_x_y, dy_omega_x_y ,ux and uy to calculate the RHS side of the equation.
Note that we also need ddx_omega_x_y and ddy_omega_x_y if viscosity is non zero.*/
	
	/*lets calculate dx_omega_x_y and dy_omega_x_y first.*/
	
	double time;
	double time_max=1;
	double time_min=0;
	double dt=0.1;
	
	#pragma omp parallel for
	for(time=time_min+dt;time<=time_max;time=time+dt)
	{  
	fftw_plan a1;
	a1 = fftw_plan_dft_r2c_2d(N1,N2,omega_x_y_o,omega_kx_ky_o,FFTW_ESTIMATE);
	fftw_execute(a1);
	fftw_destroy_plan(a1);
	
	fftw_complex *dkx_omega_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);/* derivative of phi in fourier space w.r.t kx*/
	fftw_complex *dky_omega_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);/* derivative of phi in fourier space w.r.t ky*/
	double *dx_omega_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
	double *dy_omega_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
	
	#pragma omp parallel for private(n2)
	for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  double kx = n1*dx;
			  double ky = n2*dy;
			  
			  dkx_omega_kx_ky_o[n1*N2+n2][0] = kx*(omega_kx_ky_o[n1*N2+n2][0]);
			  dkx_omega_kx_ky_o[n1*N2+n2][1] = kx*(omega_kx_ky_o[n1*N2+n2][1]);
			  /*F(d/dx(f(x))) = i*k*F(f(x)). We have to multiply by i now.We can easily swap the 
			   values of real and imaginary with a minus sign.i*(a+ib) = -b + ia.*/
			  double temp = dkx_omega_kx_ky_o[n1*N2+n2][0];
			  dkx_omega_kx_ky_o[n1*N2+n2][0] = -dkx_omega_kx_ky_o[n1*N2+n2][1]; /*real part = - imaginary part.*/
			  dkx_omega_kx_ky_o[n1*N2+n2][1] = temp;
			  
			  /* we will do the same for dky_phi_kx_ky.*/
			  dky_omega_kx_ky_o[n1*N2+n2][0] = ky*omega_kx_ky_o[n1*N2+n2][0];
			  dky_omega_kx_ky_o[n1*N2+n2][1] = ky*omega_kx_ky_o[n1*N2+n2][1];
			  double temp1 = dky_omega_kx_ky_o[n1*N2+n2][0];
			  dky_omega_kx_ky_o[n1*N2+n2][0] = -dky_omega_kx_ky_o[n1*N2+n2][1]; /*real part = - imaginary part.*/
			  dky_omega_kx_ky_o[n1*N2+n2][1] = temp1;
		  }
	  }
	  /*we will now take inverse fourier transform of these to get dx_phi_x_y and dy_phi_x_y.*/
	  fftw_plan a = fftw_plan_dft_c2r_2d(N1,N2,dkx_omega_kx_ky_o,dx_omega_x_y_o,FFTW_ESTIMATE);
	  fftw_plan b = fftw_plan_dft_c2r_2d(N1,N2,dky_omega_kx_ky_o,dy_omega_x_y_o,FFTW_ESTIMATE);
	  fftw_execute(a);
	  fftw_execute(b);
	  fftw_destroy_plan(a);
	  fftw_destroy_plan(b);
	  /*Note that we have to normalise after inverse fourier transform.*/
	  
	  #pragma omp parallel for private(n2)
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  dx_omega_x_y_o[n1*N2+n2] = dx_omega_x_y_o[n1*N2+n2]/(float(N1)*float(N2));
			  dy_omega_x_y_o[n1*N2+n2] = dy_omega_x_y_o[n1*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	/*we have calculated dx_omega_x_y_o and dy_omega_x_y_o.We are left with ddx_omega_x_y_o, 
	 ddy_omega_x_y_o,ux and uy lets calculate them.*/
	fftw_complex *ddkx_omega_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);/* derivative of phi in fourier space w.r.t kx*/
	fftw_complex *ddky_omega_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);/* derivative of phi in fourier space w.r.t ky*/
	double *ddx_omega_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
	double *ddy_omega_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
	
	#pragma omp parallel for private(n2)
	for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  double kx = n1*dx;
			  double ky = n2*dy;
			  
			  ddkx_omega_kx_ky_o[n1*N2+n2][0] = kx*dkx_omega_kx_ky_o[n1*N2+n2][0];
			  ddkx_omega_kx_ky_o[n1*N2+n2][1] = kx*dkx_omega_kx_ky_o[n1*N2+n2][1];
			  /*F(d/dx(f(x))) = i*k*F(f(x)). We have to multiply by i now.We can easily swap the 
			   values of real and imaginary with a minus sign.i*(a+ib) = -b + ia.*/
			  double temp = ddkx_omega_kx_ky_o[n1*N2+n2][0];
			  ddkx_omega_kx_ky_o[n1*N2+n2][0] = -ddkx_omega_kx_ky_o[n1*N2+n2][1]; /*real part = - imaginary part.*/
			  ddkx_omega_kx_ky_o[n1*N2+n2][1] = temp;
			  
			  /* we will do the same for dky_phi_kx_ky.*/
			  ddky_omega_kx_ky_o[n1*N2+n2][0] = ky*dky_omega_kx_ky_o[n1*N2+n2][0];
			  dky_omega_kx_ky_o[n1*N2+n2][1] = ky*dky_omega_kx_ky_o[n1*N2+n2][1];
			  double temp1 = ddky_omega_kx_ky_o[n1*N2+n2][0];
			  ddky_omega_kx_ky_o[n1*N2+n2][0] = -ddky_omega_kx_ky_o[n1*N2+n2][1]; /*real part = - imaginary part.*/
			  ddky_omega_kx_ky_o[n1*N2+n2][1] = temp1;
		  }
	  }
	  /*we will now take inverse fourier transform of these to get dx_phi_x_y and dy_phi_x_y.*/
	  fftw_plan c = fftw_plan_dft_c2r_2d(N1,N2,ddkx_omega_kx_ky_o,ddx_omega_x_y_o,FFTW_ESTIMATE);
	  fftw_plan d = fftw_plan_dft_c2r_2d(N1,N2,ddky_omega_kx_ky_o,ddy_omega_x_y_o,FFTW_ESTIMATE);
	  fftw_execute(c);
	  fftw_execute(d);
	  fftw_destroy_plan(c);
	  fftw_destroy_plan(d);
	  /*Note that we have to normalise after inverse fourier transform.*/
	  
	  #pragma omp parallel for private(n2)
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  ddx_omega_x_y_o[n1*N2+n2] = ddx_omega_x_y_o[n1*N2+n2]/(float(N1)*float(N2));
			  ddy_omega_x_y_o[n1*N2+n2] = ddy_omega_x_y_o[n1*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	/* we are now left with ux and uy calculation.Lets do them.*/
	
	
	/* we will now find phi from omega_kx_ky.*/
	/*we will find phi_kx_ky by dividing omega_kx_ky by (kx^2 + ky^2).*/
	fftw_complex *phi_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);
	
	#pragma omp parallel for private(n2)
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
			   phi_kx_ky_o[n1*N2+n2][0] = 0.0;
			   phi_kx_ky_o[n1*N2+n2][1] = 0.0; 
		   }
           else
           {
			  double temp = kx*kx + ky*ky;
			  phi_kx_ky_o[n1*N2+n2][0] = omega_kx_ky_o[n1*N2+n2][0]/temp;
			  phi_kx_ky_o[n1*N2+n2][1] = omega_kx_ky_o[n1*N2+n2][1]/temp;
		   }
		}
	}
	
/*we now have phi_kx_ky in our hand. We will do inverse fourier transform and go back to phi_x_y */

      double *phi_x_y_o;
      phi_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
      fftw_plan q;
      q = fftw_plan_dft_c2r_2d(N1,N2,phi_kx_ky_o,phi_x_y_o,FFTW_ESTIMATE);
      fftw_execute(q);
      /*Note that we have to normalise after inverse fourier transform.*/
      
      #pragma omp parallel for private(n2)
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  phi_x_y_o[n1*N2+n2] = phi_x_y_o[n1*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	  
	  /* now from this phi_x_y we have to calculate ux and uy. we need dx_phi_x_y and dy_phi_x_y 
	   for this. */
	  
	  fftw_complex *dkx_phi_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);/* derivative of phi in fourier space w.r.t kx*/
	  fftw_complex *dky_phi_kx_ky_o = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N1*N2);/* derivative of phi in fourier space w.r.t ky*/
	  double *dx_phi_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
	  double *dy_phi_x_y_o = (double*)malloc(sizeof(double)*N1*N2);
	  
	  #pragma omp for private(n2)
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  double kx = n1*dx;
			  double ky = n2*dy;
			  
			  dkx_phi_kx_ky_o[n1*N2+n2][0] = kx*phi_kx_ky_o[n1*N2+n2][0];
			  dkx_phi_kx_ky_o[n1*N2+n2][1] = kx*phi_kx_ky_o[n1*N2+n2][1];
			  /*F(d/dx(f(x))) = i*k*F(f(x)). We have to multiply by i now.We can easily swap the 
			   values of real and imaginary with a minus sign.i*(a+ib) = -b + ia.*/
			  double temp = dkx_phi_kx_ky_o[n1*N2+n2][0];
			  dkx_phi_kx_ky_o[n1*N2+n2][0] = -dkx_phi_kx_ky_o[n1*N2+n2][1]; /*real part = - imaginary part.*/
			  dkx_phi_kx_ky_o[n1*N2+n2][1] = temp;
			  
			  /* we will do the same for dky_phi_kx_ky.*/
			  dky_phi_kx_ky_o[n1*N2+n2][0] = ky*phi_kx_ky_o[n1*N2+n2][0];
			  dky_phi_kx_ky_o[n1*N2+n2][1] = ky*phi_kx_ky_o[n1*N2+n2][1];
			  double temp1 = dky_phi_kx_ky_o[n1*N2+n2][0];
			  dky_phi_kx_ky_o[n1*N2+n2][0] = -dky_phi_kx_ky_o[n1*N2+n2][1]; /*real part = - imaginary part.*/
			  dky_phi_kx_ky_o[n1*N2+n2][1] = temp1;
		  }
	  }
	  /*we will now take inverse fourier transform of these to get dx_phi_x_y and dy_phi_x_y.*/
	  fftw_plan r = fftw_plan_dft_c2r_2d(N1,N2,dkx_phi_kx_ky_o,dx_phi_x_y_o,FFTW_ESTIMATE);
	  fftw_plan s = fftw_plan_dft_c2r_2d(N1,N2,dky_phi_kx_ky_o,dy_phi_x_y_o,FFTW_ESTIMATE);
	  fftw_execute(r);
	  fftw_execute(s);
	  fftw_destroy_plan(r);
	  fftw_destroy_plan(s);
	  /*Note that we have to normalise after inverse fourier transform.*/
	  
	  #pragma omp for private(n2)
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  dx_phi_x_y_o[n1*N2+n2] = dx_phi_x_y_o[n1*N2+n2]/(float(N1)*float(N2));
			  dy_phi_x_y_o[n1*N2+n2] = dy_phi_x_y_o[n1*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	  
	  /* we have to now find ux and uy from dx_phi_x_y and dy_phi_x_y*/
	  double *ux_o = (double*)malloc(sizeof(double)*N1*N2);
	  double *uy_o = (double*)malloc(sizeof(double)*N1*N2);
	  
	  #pragma omp for private(n2)
	  for(int n1=0;n1<N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  ux_o[n1*N2+n2] = +dy_phi_x_y_o[n1*N2+n2];
			  uy_o[n1*N2+n2] = -dx_phi_x_y_o[n1*N2+n2];
		  }
	  }
	 
	 
	 /* all the values required to calculate the RHS are found.Lets call this RHS term as force
	   term.*/
	   double *force_o = (double*)malloc(sizeof(double)*N1*N2);
	   double viscosity = 1e-3;
	   
	   #pragma omp for private(n2)
	   for(int n1=0;n1<N1;n1++)
	   {
		   for(int n2=0;n2<N2;n2++)
		   {
			   force_o[n1*N1+n2] = viscosity*(ddx_omega_x_y_o[n1*N1+n2] + ddy_omega_x_y_o[n1*N1+n2]) - ux_o[n1*N1+n2]*dx_omega_x_y_o[n1*N1+n2] - uy_o[n1*N1+n2]*dy_omega_x_y_o[n1*N1+n2];
			   omega_x_y[n1*N1+n2] = omega_x_y_o[n1*N1+n2] + force_o[n1*N1+n2]*dt;
		   }
	   }
	   
	   /*lets print ux uy and time.*/
	   cout<<"time:"<<time<<endl;
	   for(int n1=0;n1<N1;n1++)
	   {
		   for(int n2=0;n2<N2;n2++)
		   {
			   cout<<n1<<" "<<n2<<" "<<"  "<<omega_x_y[n1*N2+n2]<<endl;
		   }
	   }
	   
	   /* we will now copy omega_x_y into omega_x_y_o*/
	   #pragma omp for private(n2)
	   for(int n1=0;n1<N1;n1++)
	   {
		   for(int n2 = 0;n2<N2;n2++)
		   {
			   omega_x_y_o[n1*N1+n2] = omega_x_y[n1*N1+n2];
		   }
	   }
	}
	return 0;
}
