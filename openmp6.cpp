#include<bits/stdc++.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<omp.h>
using namespace std;

int main()
{
	int sum = 0;
	int i;
	#pragma omp parallel for num_threads(2)
	for(i=0;i<=10;i++)
	{
		sum = sum + i;
		cout<<"thread id:"<<omp_get_thread_num()<<" "<<"loop:"<<i<<" "<<sum<<endl;
	}
	cout<<sum<<endl;
}
