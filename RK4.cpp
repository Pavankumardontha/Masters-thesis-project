#include<bits/stdc++.h>
#include<iostream>
using namespace std;

double f(double x,double y)
{
	double res ;
	double Num = (y*y)-(x*x);
	double deno = (y*y + x*x);
	res = Num/deno;
	return res;
}


int main()
{
	double h;
	double x0=0;
	double y0=1;
	/*we are required to find the value of y at x=0.2 and x=0.4*/
    h = 0.05;
    double n_steps = (0.4-x0)/h;
    int n = (int)n_steps;
    double x[n+1];
    double y[n+1];
    x[0]=x0;
    y[0]=y0;
    
    double k1;
    double k2;
    double k3;
    double k4;
     
    for(int i=1;i<=n;i++)
    {   
		x[i]=x[i-1]+h;
		k1 = h*f(x[i-1],y[i-1]);
		k2 = h*f(x[i-1]+h/2.0,y[i-1]+k1/2.0);
		k3 = h*f(x[i-1]+h/2.0,y[i-1]+k2/2.0);
		k4 = h*f(x[i-1]+h,y[i-1]+k3);
		y[i] = y[i-1] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
		cout<<k1<<k2<<k3<<k4<<y[i]<<endl;
	}
	cout<<"value of y at x=0.2:"<<y[4]<<endl;
	cout<<"value of y at x=0.4:"<<y[8]<<endl;
	
}
