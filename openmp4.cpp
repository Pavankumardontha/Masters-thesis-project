#include<bits/stdc++.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<omp.h>
using namespace std;

int main()
{
	/* parellising double for loops */
	int n1;
	n1 = 5;
	int n2 = 5;
	int a[5][5];
	int i;
	int j;
	#pragma omp parallel for num_threads(3) private(j)
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			cout<<"thread no; "<<omp_get_thread_num()<<" "<<i<<" "<<j<<" "<<endl;
		}
	}
	return 0;
}
