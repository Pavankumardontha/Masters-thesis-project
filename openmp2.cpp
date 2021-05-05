#include<bits/stdc++.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<omp.h>
using namespace std;

int main()
{
	int n;
	n = 5;
	int a[5];
	#pragma omp parallel num_threads(2)
	{
		for(int i=0;i<5;i++)
		{
			int id = omp_get_thread_num();
			a[i] = i;
			cout<<i<<" "<<a[i]<<" "<<id<<endl;
		}
	}
	/* The loop runs on 2 threads. The loop is not splitted. So it runs from i=0 to i=4 on both 
	   the threads. In other words it runs 10 times. How to use 2 threads and make it run only 
	   once ??
	 */
}
