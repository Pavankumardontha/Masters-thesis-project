#include<bits/stdc++.h>
#include<iostream>
#include<stdlib.h>
#include<fftw3.h>
#include<stdio.h>
#include<complex>
#include<iomanip>
#include<fftw3_mpi.h>
#include<mpi.h>
using namespace std;

int main(int argc,char* argv[])
{
	const ptrdiff_t N1,N2;
	ptrdiff_t l_N1,l_N2,l_s,n1,n2,lsize;
	/*
	  l_N1 = local N1 
	  l_N2 = local N2 
	  l_s  = local start point 
	  lsize = local size 
	  n1,n2 used for traversing our current size elements
	 */
	/*initialisation steps.*/
	MPI::Init(argc,argv);
	fftw_mpi_init();
	
	auto mpisize=MPI::COMM_WORLD.Get_size();
	auto mpirank=MPI::COMM_WORLD.Get_rank();
	
	double pi = 3.1428;
	double L = 2*pi;
	double dx = L/N1;
	double dy = L/N2;
	
	//calculating the local size
	lsize = fftw_mpi_local_size(N1,N2,MPI_COMM_WORLD,&l_N1,&l_s);
     	
    //allocating the space  	
	fftw_complex *omega_kx_ky_o;
	/* we will calculate the local space and allocate them respectively.*/
	omega_x_y_o = fftw_alloc_double(lsize);
	omega_x_y = fftw_alloc_double(lsize);
	omega_kx_ky_o = fftw_alloc_complex(lsize);
	
	/* initialising only the data corresponding to this cluster.*/
	for(n1=0;n1<l_N1;n1++)
	{
		for(n2=0;n2<N2;n2++)
		{
			omega_x_y_o[n1*N2+n2] = sin((l_s+n1)*dx)*cos(n2*dy);
		}
	}
	double time;
	double time_max=1; //0.8
	double time_min=0;
	double dt=0.1;
	
	for(time=time_min+dt;time<=time_max;time=time+dt)
	{
	fftw_plan a1;
a1 = fftw_mpi_plan_dft_r2c_2d(N1,N2,omega_x_y_o,omega_kx_ky_o,MPI_COMM_WORLD,FFTW_ESTIMATE);
	fftw_execute(a1);
	fftw_destroy_plan(a1);
	
	fftw_complex *dkx_omega_kx_ky_o = fftw_alloc_complex(lsize);
	/* derivative of phi in fourier space w.r.t kx*/
	fftw_complex *dky_omega_kx_ky_o = fftw_alloc_complex(lsize);
	/* derivative of phi in fourier space w.r.t ky*/
	double *dx_omega_x_y_o = fftw_alloc_double(lsize);
	double *dy_omega_x_y_o = fftw_alloc_double(lsize);
	for(n1=0;n1<l_N1;n1++)
	  {
		  for(n2=0;n2<N2;n2++)
		  {
			  double kx =(l_s+n1)*dx;
			  double ky = n2*dy;
			  
dkx_omega_kx_ky_o[(l_s+n1)*N2+n2][0] = kx*(omega_kx_ky_o[(l_s+n1)*N2+n2][0]);
dkx_omega_kx_ky_o[(l_s+n1)*N2+n2][1] = kx*(omega_kx_ky_o[(l_s+n1)*N2+n2][1]);
/*F(d/dx(f(x))) = i*k*F(f(x)). We have to multiply by i now.We can easily swap the 
values of real and imaginary with a minus sign.i*(a+ib) = -b + ia.*/
double temp = dkx_omega_kx_ky_o[(ls_n1)*N2+n2][0];
dkx_omega_kx_ky_o[(ls+n1)*N2+n2][0] = -dkx_omega_kx_ky_o[(l_s+n1)*N2+n2][1]; 
/*real part = - imaginary part.*/
dkx_omega_kx_ky_o[(ls+n1)n1*N2+n2][1] = temp;
			  
/* we will do the same for dky_phi_kx_ky.*/
dky_omega_kx_ky_o[(l_s+n1)*N2+n2][0] = ky*omega_kx_ky_o[(l_s+n1)*N2+n2][0];
dky_omega_kx_ky_o[(l_s+n1)*N2+n2][1] = ky*omega_kx_ky_o[(l_s+n1)*N2+n2][1];
double temp1 = dky_omega_kx_ky_o[(l_s+n1)*N2+n2][0];
dky_omega_kx_ky_o[(l_s+n1)*N2+n2][0] = -dky_omega_kx_ky_o[(l_s+n1)*N2+n2][1]; 
/*real part = - imaginary part.*/
dky_omega_kx_ky_o[(l_s+n1)*N2+n2][1] = temp1;
		  }
	  }
	/*we will now take inverse fourier transform of these to get dx_phi_x_y and dy_phi_x_y.*/
fftw_plan a = fftw_mpi_plan_dft_c2r_2d(N1,N2,dkx_omega_kx_ky_o,dx_omega_x_y_o,MPI_COMM_WORLD,
FFTW_ESTIMATE);
fftw_plan b = fftw_mpi_plan_dft_c2r_2d(N1,N2,dky_omega_kx_ky_o,dy_omega_x_y_o,MPI_COMM_WORLD,
FFTW_ESTIMATE);
	  fftw_execute(a);
	  fftw_execute(b);
	  fftw_destroy_plan(a);
	  fftw_destroy_plan(b);
	  /*Note that we have to normalise after inverse fourier transform.*/
	  for(n1=0;n1<l_N1;n1++)
	  {
		  for(n2=0;n2<N2;n2++)
		  {
	dx_omega_x_y_o[(l_s+n1)*N2+n2] = dx_omega_x_y_o[(l_s+n1)*N2+n2]/(float(N1)*float(N2));
	dy_omega_x_y_o[(l_s+n1)*N2+n2] = dy_omega_x_y_o[(l_s+n1)*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	/*we have calculated dx_omega_x_y_o and dy_omega_x_y_o.We are left with ddx_omega_x_y_o, 
	 ddy_omega_x_y_o,ux and uy lets calculate them.*/
	fftw_complex *ddkx_omega_kx_ky_o = fftw_alloc_complex(lsize);
	/* derivative of phi in fourier space w.r.t kx*/
	fftw_complex *ddky_omega_kx_ky_o = fftw_alloc_complex(lsize);
	/* derivative of phi in fourier space w.r.t ky*/
	double *ddx_omega_x_y_o = fftw_alloc_double(lsize);
	double *ddy_omega_x_y_o = fftw_alloc_double(lsize);
	for(int n1=0;n1<l_N1;n1++)
	  {
		  for(int n2=0;n2<N2;n2++)
		  {
			  double kx = (l_s+n1)*dx;
			  double ky = n2*dy;
			  
ddkx_omega_kx_ky_o[(l_s+n1)*N2+n2][0] = kx*dkx_omega_kx_ky_o[(l_s+n1)*N2+n2][0];
ddkx_omega_kx_ky_o[(l_s+n1)*N2+n2][1] = kx*dkx_omega_kx_ky_o[(l_s+n1)*N2+n2][1];
/*F(d/dx(f(x))) = i*k*F(f(x)). We have to multiply by i now.We can easily swap the 
values of real and imaginary with a minus sign.i*(a+ib) = -b + ia.*/
double temp = ddkx_omega_kx_ky_o[(l_s+n1)*N2+n2][0];
ddkx_omega_kx_ky_o[(l_s+n1)*N2+n2][0] = -ddkx_omega_kx_ky_o[(l_s+n1)*N2+n2][1]; 
			  /*real part = - imaginary part.*/
ddkx_omega_kx_ky_o[(l_s+n1)*N2+n2][1] = temp;
			  
			  /* we will do the same for dky_phi_kx_ky.*/
ddky_omega_kx_ky_o[(l_s+n1)*N2+n2][0] = ky*dky_omega_kx_ky_o[(l_s+n1)*N2+n2][0];
dky_omega_kx_ky_o[(l_s+n1)*N2+n2][1] = ky*dky_omega_kx_ky_o[(l_s+n1)*N2+n2][1];
double temp1 = ddky_omega_kx_ky_o[(l_s+n1)*N2+n2][0];
ddky_omega_kx_ky_o[(l_s+n1)*N2+n2][0] = -ddky_omega_kx_ky_o[(l_s+n1)*N2+n2][1]; 
			  /*real part = - imaginary part.*/
ddky_omega_kx_ky_o[(l_s+n1)*N2+n2][1] = temp1;
		  }
	  }
	  /*we will now take inverse fourier transform of these to get dx_phi_x_y and dy_phi_x_y.*/
fftw_plan c = fftw_mpi_plan_dft_c2r_2d(N1,N2,ddkx_omega_kx_ky_o,ddx_omega_x_y_o,MPI_COMM_WORLD,
FFTW_ESTIMATE);
fftw_plan d = fftw_mpi_plan_dft_c2r_2d(N1,N2,ddky_omega_kx_ky_o,ddy_omega_x_y_o,MPI_COMM_WORLD,
FFTW_ESTIMATE);
	  fftw_execute(c);
	  fftw_execute(d);
	  fftw_destroy_plan(c);
	  fftw_destroy_plan(d);
	  /*Note that we have to normalise after inverse fourier transform.*/
	  for(n1=0;n1<l_N1;n1++)
	  {
		  for(n2=0;n2<N2;n2++)
		  {
ddx_omega_x_y_o[(l_s+n1)*N2+n2] = ddx_omega_x_y_o[(l_s+n1)*N2+n2]/(float(N1)*float(N2));
ddy_omega_x_y_o[(l_s+n1)*N2+n2] = ddy_omega_x_y_o[(l_s+n1)*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	/* we are now left with ux and uy calculation.Lets do them.*/
	
	
	/* we will now find phi from omega_kx_ky.*/
	/*we will find phi_kx_ky by dividing omega_kx_ky by (kx^2 + ky^2).*/
	fftw_complex *phi_kx_ky_o = fftw_alloc_complex(lsize);
	
	for(n1=0;n1<l_N1;n1++)
	{
		for(n2=0;n2<N2;n2++)
		{
/* we will first calculate kx and ky. if both kx and ky equals 0 then then phi_kx_ky 
becomes infinite in this case since kx*kx + ky*ky = 0. we assume phi_kx_ky=0 when 
both kx and ky are zero.*/
		   double kx = (l_s+n1)*dx;
		   double ky = n2*dy;
           if(kx==0 && ky==0)
           {
			   phi_kx_ky_o[(l_s+n1)*N2+n2][0] = 0.0;
			   phi_kx_ky_o[(l_s+n1)*N2+n2][1] = 0.0; 
		   }
           else
           {
			  double temp = kx*kx + ky*ky;
	phi_kx_ky_o[(l_s+n1)*N2+n2][0] = omega_kx_ky_o[(l_s+n1)*N2+n2][0]/temp;
	phi_kx_ky_o[(l_s+n1)*N2+n2][1] = omega_kx_ky_o[(l_s+n1)*N2+n2][1]/temp;
		   }
		}
	}
	
/*we now have phi_kx_ky in our hand. We will do inverse fourier transform and go 
 * back to phi_x_y */

      double *phi_x_y_o;
      phi_x_y_o = fftw_alloc_double(lsize);
      fftw_plan q;
      q = fftw_mpi_plan_dft_c2r_2d(N1,N2,phi_kx_ky_o,phi_x_y_o,MPI_COMM_WORLD,FFTW_ESTIMATE);
      fftw_execute(q);
      /*Note that we have to normalise after inverse fourier transform.*/
	  for(n1=0;n1<l_N1;n1++)
	  {
		  for(n2=0;n2<N2;n2++)
		  {
			phi_x_y_o[(l_s+n1)*N2+n2] = phi_x_y_o[(l_s+n1)*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	  
	  
	  /* now from this phi_x_y we have to calculate ux and uy. we need dx_phi_x_y 
	   and dy_phi_x_y for this. */
	  
	  fftw_complex *dkx_phi_kx_ky_o = fftw_alloc_complex(lsize);
	  /* derivative of phi in fourier space w.r.t kx*/
	  fftw_complex *dky_phi_kx_ky_o = fftw_alloc_complex(lsize);
	  /* derivative of phi in fourier space w.r.t ky*/
	  double *dx_phi_x_y_o = fftw_alloc_double(lsize);
	  double *dy_phi_x_y_o = fftw_alloc_double(lsize);
	  
	  for(n1=0;n1<l_N1;n1++)
	  {
		  for(n2=0;n2<N2;n2++)
		  {
			  double kx = (l_s+n1)*dx;
			  double ky = n2*dy;
			  
			  dkx_phi_kx_ky_o[(l_s+n1)*N2+n2][0] = kx*phi_kx_ky_o[(l_s+n1)*N2+n2][0];
			  dkx_phi_kx_ky_o[(l_s+n1)*N2+n2][1] = kx*phi_kx_ky_o[(l_s+n1)*N2+n2][1];
/*F(d/dx(f(x))) = i*k*F(f(x)). We have to multiply by i now.We can easily swap the values
of real and imaginary with a minus sign.i*(a+ib) = -b + ia.*/
			  double temp = dkx_phi_kx_ky_o[(l_s+n1)*N2+n2][0];
			  dkx_phi_kx_ky_o[(l_s+n1)*N2+n2][0] = -dkx_phi_kx_ky_o[(l_s+n1)*N2+n2][1];
			   /*real part = - imaginary part.*/
			  dkx_phi_kx_ky_o[(l_s+n1)*N2+n2][1] = temp;
			  
			  /* we will do the same for dky_phi_kx_ky.*/
			  dky_phi_kx_ky_o[(l_s+n1)*N2+n2][0] = ky*phi_kx_ky_o[(l_s+n1)*N2+n2][0];
			  dky_phi_kx_ky_o[(l_s+n1)*N2+n2][1] = ky*phi_kx_ky_o[(l_s+n1)*N2+n2][1];
			  double temp1 = dky_phi_kx_ky_o[(l_s+n1)*N2+n2][0];
			  dky_phi_kx_ky_o[(l_s+n1)*N2+n2][0] = -dky_phi_kx_ky_o[(l_s+n1)*N2+n2][1]; 
			  /*real part = - imaginary part.*/
			  dky_phi_kx_ky_o[(l_s+n1)*N2+n2][1] = temp1;
		  }
	  }
	  /*we will now take inverse fourier transform of these to get dx_phi_x_y and dy_phi_x_y.*/
	  fftw_plan r = fftw_mpi_plan_dft_c2r_2d(N1,N2,dkx_phi_kx_ky_o,dx_phi_x_y_o,MPI_COMM_WORLD,
	  FFTW_ESTIMATE);
	  fftw_plan s = fftw_mpi_plan_dft_c2r_2d(N1,N2,dky_phi_kx_ky_o,dy_phi_x_y_o,MPI_COMM_WORLD,
	  FFTW_ESTIMATE);
	  fftw_execute(r);
	  fftw_execute(s);
	  fftw_destroy_plan(r);
	  fftw_destroy_plan(s);
	  /*Note that we have to normalise after inverse fourier transform.*/
	  for(n1=0;n1<l_N1;n1++)
	  {
		  for(n2=0;n2<N2;n2++)
		  {
dx_phi_x_y_o[(l_s+n1)*N2+n2] = dx_phi_x_y_o[(l_s+n1)*N2+n2]/(float(N1)*float(N2));
dy_phi_x_y_o[(l_s+n1)*N2+n2] = dy_phi_x_y_o[(l_s+n1)*N2+n2]/(float(N1)*float(N2));
		  }
	  }
	  
	  /* we have to now find ux and uy from dx_phi_x_y and dy_phi_x_y*/
	  double *ux_o = fftw_alloc_double(lsize);
	  double *uy_o = fftw_alloc_double(lsize);
	  for(n1=0;n1<l_N1;n1++)
	  {
		  for(n2=0;n2<N2;n2++)
		  {
			  ux_o[(l_s+n1)*N2+n2] = +dy_phi_x_y_o[(l_s+n1)*N2+n2];
			  uy_o[(l_s+n1)*N2+n2] = -dx_phi_x_y_o[(l_s+n1)*N2+n2];
		  }
	  }
	 
	 
	 /* all the values required to calculate the RHS are found.Lets call this RHS term as 
	  force term.*/
	   double *force_o = fftw_alloc_double(lsize);
	   double viscosity = 1e-3;
	   for(n1=0;n1<l_N1;n1++)
	   {
		   for(n2=0;n2<N2;n2++)
		   {
force_o[(l_s+n1)*N1+n2] = viscosity*(ddx_omega_x_y_o[(l_s+n1)*N1+n2] + 
ddy_omega_x_y_o[(l_s+n1)*N1+n2]) - ux_o[(l_s+n1)*N1+n2]*dx_omega_x_y_o[(l_s+n1)*N1+n2] 
- uy_o[(l_s+n1)*N1+n2]*dy_omega_x_y_o[(l_s+n1)*N1+n2];
omega_x_y[(l_s+n1)*N1+n2] = omega_x_y_o[(l_s+n1)*N1+n2] + force_o[(l_s+n1)*N1+n2]*dt;
		   }
	   }
	   
	   /*lets print ux uy and time.*/
	   cout<<"time:"<<time<<endl;
	   cout<<"n1"<<"\t"<<"n2"<<"\t"<<"omega_z"<<endl;
	   for(n1=0;n1<l_N1;n1++)
	   {
		   for(n2=0;n2<N2;n2++)
		   {
			   cout<<(l_s+n1)<<"\t"<<n2<<"\t"<<omega_x_y[(l_s+n1)*N2+n2]<<endl;
		   }
	   }
	   
	   /* we will now copy omega_x_y into omega_x_y_o*/
	   for(n1=0;n1<l_N1;n1++)
	   {
		   for(n2 = 0;n2<N2;n2++)
		   {
			   omega_x_y_o[(l_s+n1)*N1+n2] = omega_x_y[(l_s+n1)*N1+n2];
		   }
	   }
	}
	MPI::finalize();
	return 0;
}
