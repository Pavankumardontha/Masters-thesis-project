#include<bits/stdc++.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<omp.h>
using namespace std;

int main()
{
	/* parellising single for loop */
	int n;
	n = 10;
	int a[10];
	int i;
	#pragma omp parallel for num_threads(3)
	for(i=0;i<10;i++)
	{
			int id = omp_get_thread_num();
			a[i] = i;
			cout<<"Thread id:"<<id<<" "<<"loop index:"<<i<<" "<<a[i]<<" "<<endl;
	}
	
}
