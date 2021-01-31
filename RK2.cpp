#include<bits/stdc++.h>
#include<iostream>
using namespace std;
/*
 given 
 1)dy/dx as a function of (x,y) ie dy/dx = f(x,y) and
 2)initial value of y,ie y(0)=some value 
 Our task is to find the unknown function y at a given point x ie y(x). 
  
  ex1: dy/dx = 3x^(2)+y^(2)
  find the value of y at x = 1.1 given y=1.2 @ x=1 
  
  solution: 
  RK2 formulae:-
  k1 = hf(x[n],y[n]) 
  k2 = hf(x[n]+h,y[n]+k1)
  y[n+1]=y[n]+1/2(k1+k2)
  
  given x0=1 and y0=1.2 h=0.1 
  k1=hf(x0,y0)=0.1(4.44)=0.444
  k2=hf(x0+h,y0+k1)=0.1*f(1+0.1,1.2+0.444) = 0.600
  thus k1 = 0.444 and k2=0.600
  y1=y0+1/2(k1+k2) = 1.2 + 1/2(0.444+0.600)= 1.722
  */
   
double f(double x,double y)
{
	double res = 3*x*x + y*y;
	return res;
}
  
  
 
int main()
{
	double h=0.05;
	double x0=1.00;
	double y0=1.2;
	double n_steps = (xn-x0)/h;
	int n = (int)n_steps;
	double x[n+1];
	double y[n+1];
	x[0]=x0;
	y[0]=y0;
	double k1,k2;
	for(int i=1;i<=n;i++)
	{
		x[i]=x[i-1]+h;
		k1=h*f(x[i-1],y[i-1]);
		k2=h*f(x[i-1]+h,y[i-1]+k1);
		y[i]=y[i-1]+0.5*(k1+k2);
	}
	cout<<"value of y at x=1.1 is:"<<y[n]<<endl;
}
