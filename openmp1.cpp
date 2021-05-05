#include<bits/stdc++.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<omp.h>
using namespace std;

int main()
{
	#pragma omp parallel num_threads(10)
	{
		int id = omp_get_thread_num();
		cout<<"Hello world from: "<<id<<endl;
	}
	return 0;
	/*run this and next one to understand the meaning of critical section */
}

