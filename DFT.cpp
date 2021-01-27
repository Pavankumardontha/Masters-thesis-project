#include<bits/stdc++.h>
#include<iostream>
#include<fftw3.h>
using namespace std;

int main()
{
   int n=16;
   fftw_complex in[n],out[n],in2[n];
   fftw_plan p,q;
   int i;
   
   /* cosine wave of frequency 3*/
   for(int i=0;i<n;i++)
   {
	   //cos(2 * pie * f * t)
	   /* if we fix the frequency,then our entire wave is fixed.we have to just vary the time and calculate the function value at these varying time values.*/
	   in[i][0] = cos(3*2*3.14*i/n); //we have a frequency of 3 here.
	   in[i][1]=0;
   }
   
   /* forward fourier tranform.saving in out */
   p = fftw_plan_dft_1d(n,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
   fftw_execute(p);
   cout<<"freq:"<<" "<<"cosine"<<" "<<"sine"<<endl;
   for(int i=0;i<n;i++)
   cout<<i<<"    "<<out[i][0]<<"    "<<out[i][1]<<endl;
   fftw_destroy_plan (p);
   fftw_cleanup ();
   return 0;
   /*
     we have to compile the program using command : 
     g++ DFT.cpp -lfftw3 -lm 
    */
}
