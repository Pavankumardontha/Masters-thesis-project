#include<bits/stdc++.h>
#include "mpi.h"
using namespace std;

/*we will write a function that we will integrate .y=x^(2) within the given bounds.*/
double func(double x)
{
	return x*x;
}

int main(int argc,char* argv[])
{
	/*the above 2 inputs will come from command line arguments.argc=no of inputs from the command line.argv[] is an input array.*/
	/*lets do initialisation.*/
	MPI::Init(argc,argv);
	
	auto mpisize=MPI::COMM_WORLD.Get_size();
	auto mpirank=MPI::COMM_WORLD.Get_rank();
	
	auto mpiroot=0;//we can take any ID as our root.Here we are considering node 0 as our root node.
	int nRectangles;
	double xLeft;
	double xRight;
	double domainWidth;
	double rectangleWidth;
	char rule;//this is the rule that we use to do integration.We will use rieman sum and tripezoidal in the next step
	
	/*lets write some code that uses the broadcast from MPI.Note that we will take all the input in only root rank.After the root rank takes
	these above mentioned inputs,it will use the broadcast function of the MPI and broadcast these values into the different roots.*/
	
	if(mpirank==mpiroot)
	{
		/*we have to take inputs in this rank of all the above mentioned variables.*/
		if(argc==5)//since we are expecting 5 input varibles from the command line arguments.if they are not 5 then we will print the error.
		{
			/*note that we have an array of char pointers given as input by the user.We have to now change these.We have to convert these
			inputs into the above datatypes.*/
			nRectangles = atoi(argv[1]);//note that we have to start from 1.This is very important.argv[0] stores the name of our file name
			xLeft=atof(argv[2]);//the left most x-coordinate
			xRight=atof(argv[3]);//the right most x-coordinate
			rule=*argv[4];//argv[4] stores the address.
			if(rule!='R' && rule!='T' && rule!='r' && rule!='t' )
			{
				cerr<<"incorrect command line argument.Must be 'R' or 'r' for rectangle method and 'T' or 't' for trapezoid "<<endl;
				exit(1);
			}
			domainWidth=xRight-xLeft;
			rectangleWidth=domainWidth/nRectangles;
/* all these values are colected on rank 0.Those values are only assigned to rank 0.On everyother rank it doesnt have values yet.
On other IDs,all these variables are junk values.So we have to broadcat these values from ID=0.Note that we have to not write
the broadcast function in this if condition.If we do so ,then only the rank 0 will broadcast the values.but in our tree model,
every ID that has got values from the root must also act as helper function and thus must broadcast the values to other IDs.
so we will broadcast the values outside this If condition so every ID executes it and not only ID 0*/		
		}
		else
		{
			cerr<<"incorrect no of command-line arguments!!!"<<endl;
			exit(1);
		}
		
	}
	  /* we have taken 4 inputs(excluding the file name in arg[0]) in the root node.From those 4 inputs,we have caculated 
	  the rectangle width and domain width so we are indirectly
	  broadcasting 6 variables to all other nodes from the root node.*/
		/*Broadcasting the values.Note that we can broadcast from the root as well as the helper functions.*/
		MPI::COMM_WORLD.Bcast(nRectangles,1,MPI::INT,mpiroot);
		/* lets do the same for everyother parameters.*/
		MPI::COMM_WORLD.Bcast(&xLeft,1,MPI::DOUBLE,mpiroot);
		MPI::COMM_WORLD.Bcast(&xRight,1,MPI::DOUBLE,mpiroot);
		MPI::COMM_WORLD.Bcast(&domainWidth,1,MPI::DOUBLE,mpiroot);//domain width is calculated based on 4 inputs
		MPI::COMM_WORLD.Bcast(&rectangleWidth,1,MPI::DOUBLE,mpiroot);//rectangle width is calculate based on 4 inputs
		MPI::COMM_WORLD.Bcast(&rule,1,MPI::CHAR,mpiroot);
		
		
		/*DISTRIBUTION OF WORK AMONG DIFFERENT NODES .*/
	/*we have to distribute the work equally so let us first calculate how many rectangles area each node has to compute.*/
		/*we are done with our broadcasting.We have to distribute the no of rectangles on each ID.*/
		auto nRectanglesLocal=nRectangles/mpisize;//this stores the no of rectangles in current ID.mpisize is the total no of ranks or IDs.
		/*some of the rectangles will miss out if we do this above method of division.We will distribute the left over rectangles into the
		last ID.*/
		if(mpirank==mpisize-1)//ranks start indexing from 0.
		{
			/*this is for distributing the remaining rectangles to the last node.*/
			nRectanglesLocal += nRectangles % mpisize;
			/*if there are 11 rectangles and mpisize is 5,then ID=0 gets 2 rectangles,ID=1 gets 2,ID=2 gets 2,ID=3 gets 2, but ID=4 gets 3
			rectangles.*/
			
		}
	
	    /*we will have different left boundary for each of the rectangles.*/
		xLeftLocal=(xLeft + rectangleWidth*mpirank*(nRectangles/mpisize);//this stores the left most boundary of the rectangles distributed on each of the ID.
		//for rank 0 its xleft;
	    //for rank 1,its the right coordinate of the starting rectangle and so on
	    
	    /* so we have computed the left boundary of starting rectangle in each ID.So lets just print something to check what we are doing is 
	    fine.After that we will comment out this.*/
	    cout<<"rank: "<<mpirank<<"nRectanglesLocal : "<<nRectanglesLocal<<" xleft Local: "<<xLeftLocal<<endl;
	    
	    auto sumLocal=0.0d;//double precision number
	    if(rule=='R' || rule== 'r')
	    {
	    	//write code for rectangles.
	    	for(int iRectangle=0;iRectangle < nRectanglesLocal;++iRectangle)
	    	{
	    		auto myleft = xleftLocal + iRectangle * rectangleWidth;
	    		auto myHeight = func(myLeft); //y-value corresponding to myLeft x-value
	    		sumLocal +=rectangleWidth * myHeight;
	    	}
	    	/*so now we have the local sum.We have to collect all those localsums and store in the global sum variable.For this we use mpi 
	    	reduction operator.*/
	    	
		}
		else if(rule=='T' || rule=='t')
		{
			//write code for computing integral using trapezoidal rule.
			//we need a vector to store the heights of the trapezoid.Basically these are the end points of the rectangle heights
			vector<double> heights(nRectanglesLocal+1,0.0);
			//this vector has all the heights embedded in it.
			for(auto iRectangle=0;iRectangle < nRectanglesLocal + 1; ++iRectangle)
			{
				auto myLeft = xLeftLocal + iRectangle * rectangleWidth; 
				heights[iRectangle] = func(myLeft); //our function is y = x^(2)
			}
			//we have our heights vector ready with us. !!!Lets apply trapezoidal formula 
			
			for(auto iRectangle=0;iRectangle < nRectanglesLocal; ++iRectangle)
			{
				auto myArea = (heights[iRectangle] + heights[iRectangle + 1])/ 2.0 * rectangleWidth;
				sumLocal += myArea;
			}
		}
	    /*Using MPI reduce to collect the local sums and sending it to the root node.*/
	    auto sum=0.0d; //we will store all the sumLocal varaibles sum into 'sum' varaible.
	    MPI::COMM_WORLD.Reduce(&sumLocal,&sum,1,MPI::DOUBLE,MPI::SUM,mpiroot); 
	    //since we are passing the mpiroot,this global sum has the result only on the root node
	    /*this global sum has the correct result only on the root.*/
	    if(mpirank == mpiroot)
	    {
	    	//since our function is x^2 we can calculate the integral of this function manually.Lets calculate and try to compare our results.
	    	auto exact=(xRight * xRight * xRight - xLeft * xLeft * xLeft)/3.0;
	    	cout<<"Aread under the curve :"<<sum<<" "<<"exact value: "<<exact<<endl
	    	exit(1);
		}
	MPI::finalize();
	return 0;
}
