#include<bits/stdc++.h>
#include<iostream>
using namespace std;
#include<fftw3.h>
#include<stdio.h>
#include<stdlib.h>

int main()
{
	#pragma omp parallel 
	{
		cout<<"hello world"<<endl;
	}
	return 0;
}

/* note that we havent included any of the openmp library headers here.
gcc implements openMP standard. We have to enable it by using command:

if we have C program:
gcc -fopenmp filename.c 

if we have cpp program:
g++ -fopenmp filename.cpp

here since our filename is openmpimport.cpp we use the command to compile:
g++ -fopenmp openmpimport.cpp  

after the compilation is just simple to execute the file using command 
./a.out


*/
