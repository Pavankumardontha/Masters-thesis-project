#include<bits/stdc++.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<omp.h>
using namespace std;

int main()
{
	int n1 = 3;
	int n2 = 3;
	int n3 = 3;
	int i;
	int j;
	int k;
	#pragma omp parallel for num_threads(3) private(j) private(k)
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			for(k=0;k<n3;k++)
			{
				cout<<"thread no:"<<omp_get_thread_num()<<" "<<i<<" "<<j<<" "<<k<<endl;
			}
		}
	}
	return 0;
}
