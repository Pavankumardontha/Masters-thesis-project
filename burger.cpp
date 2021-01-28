#include<bits/stdc++.h>
#include<iostream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<complex>
#include<fftw3.h>
#include<math.h>
//#include<rfftw.h>
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
	int N=10;
	double pi = 3.1428;
	double L;//length of the interval 
	double dx;
	double x[N];
	double *ux;//this is our function which is sin(x)
	L=2*pi; //length of the interval is 2pie
	dx = L/N; 
	fftw_plan p;
	fftw_complex *uk;
	ux = (double*)malloc(sizeof(double)*N);
	uk = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	for(int i=0;i<N;i++)
	{
		x[i] = (i)*dx;
		ux[i] = sin(x[i]);
	}
	/* step-1:lets calculate the fourier transform of our function ux.*/
	p = fftw_plan_dft_r2c_1d(N,ux,uk,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	
	/*
	 step-2:
	 we will calculate the the derivative of ux in the fourier space i.e the fourier transform 
	 of the derivative of ux.So we are taking the derivative of uk.*/
	 fftw_complex duk[N];
	 fftw_complex dduk[N];
	 complex<double> j(0,1); //mathematical i.
	 	 
	 /* duk = i*omega*uk  omega = ((2*pie)/N)*k. k varies from 0 to N-1. So we can directly 
	  calculate the value of duk from uk by juts looping through.*/
	  double W_real;
	  complex<double> temp;
	  
	  double time;
	  double time_max=1.0;
	  double time_min=0;
	  double dt=0.1;
	  
for(time=time_min;time<=time_max;time=time+dt)
   {
	  for(int i=0;i<N;i++)
	  {
		  W_real = ((2*pi)/N)*i;
		  complex<double> w(W_real,0);
		  temp=mul(j,w);
		  temp=mul(temp,complex<double>(uk[i][0],uk[i][1]));//j*w*uk
		  duk[i][0]=temp.real()/N; //division by N is for normalisation.
		  duk[i][1]=temp.imag()/N;
		  
		  /*we also need double derivative while solving the navier stokes equation.*/ 
		  temp = mul(w,w); //w*w
		  temp= -mul(temp,complex<double>(duk[i][0],duk[i][1])); //-w*w*duk
		  dduk[i][0] = temp.real()/N;
		  dduk[i][1] = temp.imag()/N;
	  }
	 
	 /* step-3: from duk[] and dduk[](derivatives in k-space) we will take the inverse fourier 
	  transforms and find the values of dux[] amd ddux[](derivatives in x-space). */
	  double* dux;
	  double* ddux;
	  fftw_plan r,q;
	  dux = (double*)malloc(sizeof(double)*N);
	  ddux = (double*)malloc(sizeof(double)*N);
	  /* note that while doing inverse fourier transform we should do normalisation as well.*/
	  r = fftw_plan_dft_c2r_1d(N,duk,dux,FFTW_ESTIMATE);
	  q = fftw_plan_dft_c2r_1d(N,dduk,ddux,FFTW_ESTIMATE);
	  fftw_execute(r);
	  fftw_execute(q);
	  fftw_destroy_plan(r);
	  fftw_destroy_plan(q);
	  fftw_cleanup();
	  
	  /* we will check once by printing our values.*/
	  
	  /*
	  cout<<"dux"<<"\t"<<setw(12)<<"ddux"<<"\t"<<endl;
	  for(int i=0;i<N;i++)
	  {
		cout<<dux[i]<<"\t"<<setw(12)<<ddux[i]<<"\t"<<endl;  
      }
      */
      
      /* we will now proceed to solve the equation since we have all the tools with us now.*/
      
      /*step-4 : defining the force term.*/
      double force[N]; //this is used to change the function ux[i] as a function of time.
      double nu = 0;
      for(int i=0;i<N;i++)
      {
		  //force[i] = -ux[i]*dux[i] + nu*ddux[i] ; nu is the viscocity coefficient
		  force[i] = -ux[i]*dux[i] + nu*ddux[i]; //note that the equation wont proceed unless and untill the value of nu(viscosity) is choosen.
		  //now the above force will update the value of ux[i] for the subsequent time.
		  ux[i] = ux[i] + force[i]*dt; 
	  }
	/* we can plot u(i,t) vs i at each step to clearly visualise how u(i,t) is evolving in time.*/  
	  
  }
	  /* 
	   so lets start to understand from the first what we have done.
* Initially we are given ux(i) at t=0 and our goal is to simulate it for other times as well. In
other words we have to find ux(i,t). 
* We have the equation 
  dut(i,t) + ux(i,t)*dux(i,t) = nu*ddux(i,t)  
=> dut(i,t) = -ux(i,t)*dux(i,t) + nu*ddux(i,t) 
let force(i,t) = force at position x and at time t = -ux(i,t)*dux(i,t) + nu*ddux(i,t)
=>dut = force 
so to update the value of u(i,t) as time progresses,we need the force at this moment of time
so we have to calculate the force term to get all the values of ux(i,t+delta(t)).
*We have to calculate 2 things to calculate the force,one is dux(i,t) and ddux(i,t).To calculate
both of these terms,we will first calculate the fourier transform of ux(i,t) which is uk(i,t).
From this we will calculate duk(i,t) and dduk(i,t). If these 2 things are known,then our job is
done almost.We will do inverse fourier transform and calculate dux(i,t) and ddux(i,t). Thus force
calculate is done. 
*/ 
	  return 0;
	
}
