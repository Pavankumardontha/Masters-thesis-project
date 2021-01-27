#include<bits/stdc++.h>
#include<iostream>
#include<math.h>
#include<complex>
#include<iomanip>
#include<stdlib.h>
#include<fftw3.h> //this is not needed here since we are writing base code for DFT
using namespace std;


vector<complex<double>> dft(vector<complex<double>> &X)
{
	//determine the no. of datapoints 
	int N = X.size();
	//if the size of the datapoints is N,then we will have N fourier transformed datapoints 
	int K = N; 
	//K is size of fourier transformed data 
	
	//we will now allocate memory for the output since we know its size
	vector<complex<double>> output;
	
	/* for calculating each fourier transform , we need to traverse whole of the input data points.
	 This is can be clearly seen from the formula which is used to obtain the fourier transformed 
	 data point from the input data points.*/ 
	
	/* for calulating the entire sum of real and imaginary we require a temporary variable.Note that we have N terms in the 
	 calculation of the fourier transformed data point.Each of the N terms is complex and we have to add this result to get 
	 the final transformed data point.*/
	 complex<double> intsum;
	
	//looping through each value of K 
	for(int k=0;k<K;k++)
	{
		intsum = complex<double>(0,0); //initialising the real and imaginary part to 0.
		//looping through each value of N
		for(int n=0;n<N;n++)
		{
		   /* we have X[n] with us we have to just calculate the exponential term and multiply both.*/
		   double realPart = cos(((2*M_PI)/N) * k * n);
		   double imagPart = sin(((2*M_PI)/N) * k * n);
		   complex<double> res(realPart,-imagPart);
		   /* we have to now multiply this term with X[n] and then add the final result to intsum.*/
		   intsum += X[n]*res;
		}
		//we have to push intsum into our output 
		output.push_back(intsum);
	}
	return output;
	
}



int main()
{
	/* we will create a test signal and calculate dft using above function.*/
	int N=1000;
	vector<complex<double>> signal;
	/* we have to fix the frequency and phase of above signal.*/
	double sigK = 3.0;
	double sigPhase = 0.0;
	
	/* we will calculate the function value at each n and store this in the input vector.*/
	
	/* we are using the function 
	  cos(((2*pie)/N)*3*n); n = [0,...N-1]
	*/
	for(int n=0;n<N;n++)
	{
	    double realPart = cos(((2*M_PI)/static_cast<double>(N))*sigK*static_cast<double>(n) + sigPhase);
	    double imagPart = 0.0;
	    complex<double> x(realPart,imagPart);
	    signal.push_back(x);
	}
	
	/*we will calculate the fourier transformed data of this input data.*/
	vector<complex<double>> output = dft(signal);
	
	/*we will print the first 6 output fourier transformed data points.*/
	cout<<"**********"<<endl;
	cout<<"first 6 samples"<<endl;
	cout<<"\n"<<"k\t"<<setw(12)<<"Real\t"<<setw(12)<<"Imag"<<endl;
	for(int i=0;i<6;i++)
	{
		cout<<i<<"\t"<<setw(12)<<output[i].real()<<"\t"<<setw(12)<<output[i].imag()<<endl;
	}
	
	return 0;
	
	
}
