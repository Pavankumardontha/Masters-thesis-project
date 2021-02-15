#include<iostream>
#include<bits/stdc++.h>
#include<omp.h>
using namespace std;

/*first study the meaning of threads,cores and processors and parallel programming.*/
/*program-1: basic hell world program in OPEN MP
int main()
{
	#pragma omp parallel
	{
		int ID=omp_get_thread_num();
		cout<<"hello world"<<" "<<ID<<endl;
	}
}

*/


/*each time we run it,the hello world comes in different order.This is because the threads are interleaved.Each time we run 
we get a different interleaving of the threads.Its the programmer job to make sure that every single way you interleave those 
threads gives the correct result.In other words,the interleaving of threads which change each time should not impact our final
result.
here the ID stores current thread id.We asked the system to give default no of threads.So the system may give us any no of 
threads and so the value of ID will not be the same everytime we run the program.
The kind of hardware the open MP assumes are Shared Memory computers. Any computer composed of multiple processing elements 
that share a common address space is called shared memory computer.These are of 2 types of multi processors that these computers
contain. So now we have multi-processors but these multi-processors share a commom memory space 
a)Symmetric Multiprocessor(SMP): A shared address space with "equal-time" access for each processor,and the OS treats every 
processor the same.

b)Non-uniform Address space Multiprocessor(NUMA): Different memory regions have different access costs.Think of memory
segmented into "near" and "far" from a processor.The memory which is "near" to a particular processor has less time access than
the memory which is far.

Modern processors such as i7 are often called an SMP.But is it really SMP? It has 6 cores and they all share a last level 
cache. But is it  true that the entire memory is equally accessible by all the processors ? The breakdown between processors
and memory  because of that cache architechture if a block you are looking at is sitting in a cache region closer to a processor,
then you have much faster acess time than if its far away.So even though we pretend that modern CPU's are SMP's they are NOT.If
something is sitting in a cache close to a processor then it has much faster access time than if its sitting in deep RAM far
away.There is no relation between threads and cores.You may have more threads than processor cores.Nothing restricts us from
having only one thread per core.All the threads run simultaneously.Its not that thread 0 ends and then only thread 1 starts 
or anything like that.Threads will interleave in a complicated way.Interleaving should be thought of as a combination.Every 
single way you interleave those threads must give you a correct result.This is the job of the programmer.This is the challenge
in multi-threaded programming. So open Mp is a multi-threaded programming language for shared address space machines.

Threads communicate by sharing variables.The data in the heap is accessible by all the threads.You got variables sitting in the heap 
of a process which can been seen and accessed by all the threads.So the threads interact by looking and accessing the variables
in the heap.Sometimes,you have data that they are sharing and this sharing doesnt occur in a SAFE way. You may have 1 thread 
thats dumping a particular variable which the other thread is trying  to access or read.The result you get changes as you 
interleave the threads in different order.We call these as RACE conditions.Unintended sharing of data causes race conditions.
Generally a condition where the program's output changes as the threads are scheduled differently each time you run the 
program is a race condtion. We control race conditions by organising and controlling access to shared variables.You do this 
with schronisation.We use shncronisation techniques to protect and order accesses to shared variables.Now 
synchronisation is really expensive.If you synchronise alot,you will get no performance.So as we write our problems our 
objective should be to synchronise as little as possible.This can be done by paying close attention to our data environment.
Open MP is used to manage the data and also how it is shared with the different threads.


NOTE ON SYNCHRONISATION:-
synchronization refers to one of two distinct but related concepts: synchronization of processes, and synchronization of data.
Process synchronization refers to the idea that multiple processes are to join up or handshake at a certain point, in order
to reach an agreement or commit to a certain sequence of action.The need for synchronization does not arise merely in 
multi-processor systems but for any kind of concurrent processes; even in single processor systems. Mentioned below are 
some of the main needs of synchronisation:

Forks and Joins: When a job arrives at a fork point, it is split into N sub-jobs which are then serviced by n tasks. After 
being serviced,each sub-job waits until all other sub-jobs are done processing. Then, they are joined again and leave the 
system. Thus, parallel programming requires synchronization as all the parallel processes wait for several other processes 
to occur.

Thread synchronization is defined as a mechanism which ensures that two or more concurrent processes or threads do not 
simultaneously execute some particular program segment known as critical section. Processes' access to critical section is 
controlled by using synchronization techniques. When one thread starts executing the critical section (serialized segment of 
the program) the other thread should wait until the first thread finishes. If proper synchronization techniques[1] are not 
applied, it may cause a race condition where the values of variables may be unpredictable and vary depending on the timings 
of context switches of the processes or threads.

For example, suppose that there are three processes, namely 1, 2, and 3. All three of them are concurrently executing, and they
need to share a common resource (critical section).Synchronization should be used here to avoid any conflicts for accessing 
this shared resource.Hence, when Process 1 and 2 both try to access that resource, it should be assigned to only one process 
at a time. If it is assigned to Process 1, the other process (Process 2) needs to wait until Process 1 frees that resource.
Another synchronization requirement which needs to be considered is the order in which particular processes or threads should
be executed.
Other than mutual exclusion, synchronization also deals with the following:
1)DEADLOCK, which occurs when many processes are waiting for a shared resource (critical section) which is being held by some 
other process. In this case, the processes just keep waiting and execute no further.
2)STARVATION, which occurs when a process is waiting to enter the critical section, but other processes monopolize the 
critical section, and the first process is forced to wait indefinitely.
3)PRIORITY INVERSION, which occurs when a high-priority process is in the critical section, and it is interrupted by a 
medium-priority process. This violation of priority rules can happen under certain circumstances and may lead to serious 
consequences in real-time systems.
4)BUSY WAITING, which occurs when a process frequently polls to determine if it has access to a critical section. This 
frequent polling robs processing time from other processes.


How do we create threads in multi-threading programming ?
The fundamental model behind OPEN MP is the fork-join Parallelism.In this the program starts as a single thread(master 
thread).It encounters with some part of the program where addition of threads coudl help out.So what it does is,it forks
a number of threads at that point.So now we have just entered from the serial part of program into parallel region.Note 
that apart from the block of threads created we also have the master thread running along.The master thread as id=0.
The other threads id no. start by numbering from 2.These collection of threads are a team of threads.Now these team of
threads work in parallel.After they are finished they have to now join back together to form a single thread again.
So after the completion we have only one thread(same like the master thread before).So now,the same one thread keeps 
running all along untill it sees another part of the program where running multiple threads makes sense.So now it again 
splits into multiple threads,each thread do its job and after the job all the threads join into a single thread.So this
is the process.You move along,fork threads,join threads and you go on doing the same.You fork while entering the parallel
regions of the program and you join when entering the sequential part of the program.In parallel region there are multiple
threads and in sequential part of the program only a single thread exists.We can we nest threads.So inside a parrallel
region,you can have a thread calling multiple other threads and so on.So nested threads are also possible in the parallel 
regime.The only way to create threads in open MP is with parallel construct.If this is not present then there are no 
multiple threads in your program.
omp_set_num_threads(4) : this command sets 4 threads

Each thread will run the code inside the pragma omp parallel.If we allocate the data outside the pragma omp parallel,this 
data sits on the heap and is visible and accessible by all the threads.If we allocate the ddata inside the pragma omp 
parallel,its allocated on the threads individual stack and is local to the allocated thread and private to all the other
threads. 

So we have to pay close attention with private and shared variables.Now lets try out some integration problem. We will 
write normal program to compute the value of the integral and then induce the parallel regime into it and see if its 
executing time is reduced. We will go through some more basic functions before doing this task

int omp_get_num_threads(): this asks for no of threads from the user

int omp_get_num_thread_num(): this gives the thread ID or rank of the thread in the parallel regime.

double omp_get_wtime() : this fucntion is used to measure the time of execution.Lets say we declare 
this before and after the pragma omp block of code,then the difference in time is the amount of time
it takes to run our parallel block of code.This lets you track the time.
Lets get into our program of computing the integral.

*/

/*program-2: to find the integral value of 4/(1+x^(2)) in the limits 0 to 1

int main()
{
	int num_steps=100000; //no of blocks
	double step=(1.0)/num_steps;//width of each block
	double x,pi,sum;
	//pi represents our final answer
	sum=0.0;//this sums all the vertical strips
	for(int i=0;i<num_steps;i++)
	{
		x=(i+0.5)*step;
		sum=sum+(4.0)/(1+x*x);
	}
	pi=sum*step;
	cout<<pi<<endl;
}

OUTPUT:-
3.14159

We will now use SPMD algorithmic strategy to compute the integral using parallel computing.This is a commonly used algorithm.
You have single program and you run it multiple times.So each thread has its own copy.Note that parallel program with one 
thread runs slower than the same program run sequentially.Parallel program with one thread is same as sequentuial program.
It runs slower because you need to form new extra variables in the parallel program.So parallel programs with one thread run
slightly slower than the sequential ones.

Each thread has its own copy but its going to use the id of the thread and get knownledge of how many total threads there are.
So each thread,from its id,knows how many total threads they are.Depending on the total no. of threads,each thread does its 
own job.So the job of a thread changes as the number of threads in the program changes.If you see the sequential prgram of 
the integral computation it has a loop which runs over 0 to 100000.By using the threads,we can split this loop.Each thread 
computes a particular part of the total loop.Note that the thread id numbering starts from 0 and goes to n-1 where n is the 
number of threads.So we use an array sum[n],instead of only one sum variable in the sequential programe to compute the work.
sum[0] represents the sum computed by thread 0.So sum[i] represents the sum calculated by thread i.We have to combine all 
these sum array elements to compute the total sum and this is our answer.We have to also split the loop work over to all the
threads such that any loop is not computed twice.Each loop must be traversed only once and in only in one thread.Example lets
say we have 4 threads with us.Then thread number 0 gets iteration 0,4,8,12,16......Thread 1 gets iteration 1,5,9,13,17,21,25.
......Thread 2 gets iterations 2,6,10,14,18,22........Thread 3 gets 3,7,11,15,19....and in this way by using 4 threads we have 
split the loop and each thread gets a certain number of iterations.Note that no iteration is computed twice in any thread.
Also a particular loop iteration is not computed in 2 threads.Each iteration is only computed once and in only any ONE thread.
We have to distribute the work following these rules.So this is cyclic or round robin policy.Each thread gets a subset of 
iterations to compute the total sum.

*/

/*program-3:-using parallel computing to compute the integral value 

int main()
{
	int num_steps=100000;
	double step=1.0/num_steps;
	int num_threads=2;//this can also be taken as input from the user.This represents the no of threads
	omp_set_num_threads(num_threads);
	double sum[num_threads];
	for(int i=0;i<num_threads;i++)
	sum[i]=0.0;
	double pi=0.0;//our final answer
	/*now these all variables are formed on the heap and can be used and modified by all the threads.The below code runs
	as many times as the no. of threads we have in this program.We have 2 here so the parallel block is executed 2 times.
	it basically defines the work of each thread.Whatever variables we define outside the pragma omp parallel is accessible
	to all the threads.These are like the common variables.They can be modified and accessed by all the threads in the 
	program.But the varibles inside the pragma omp parallel block are only accesible by the current thread and are private
	to all other threads.
	#pragma omp parallel
	{
		double x;
		int id=omp_get_thread_num(); //current thread id number
		//note that these are the 2 local variables and are on the stack.
		for(int i=id;i<num_steps;i=i+num_threads)
		{
			x=(i+0.5)*step;
			sum[id]=sum[id]+4.0/(1+x*x);
		}
	}
	/*we have to now sum up the total area computed by thread 1  and thread 2.note that the work is distributed evenly here.
	Each thread loops of around 50000 which is half of 100000.So work is distrubted almost evenly on all the threads but 
	this might not always be the case.we cannot always distribute the work equally on all the threads.Some threads do slightly
	more and some do slightly less.*/


	for(int i=0;i<num_threads;i++)
	pi=pi+sum[i]*step;
	cout<<pi<<endl;
}

/*
so at the end you have an array sum of size equal to the no of threads.Each thread provides an amount sum[i] to the total sum.
we can drastically see that by changing the no of threads the execution time decreases.There are few important points to be 
noted here.The master thread,the thread with id=0 saves the copy of the no of threads in the whole program.So using this also
you can get how many threads are there in our program.So we use

id=omp_get_thread_num();
nthrds=omp_get_num_threads();
if(id==0)
nthreads=nthrds;

note here that nthrds is a local variable and nthreads is a heap variable and should declared outside the pragma segment of the
code.So using the above method also we can find the number of threads in our program.Notice that only the master id can track
the exact number and if the id is not of the master thread then it cannot if used to find the no of threads running in our 
program.So the if condition is very important in this case.Only the master thread can find the total no of threads in our 
program.But why do we need to find out the no of threads using master thread ??
We you enter parallel region in open MP you request a no of threads.But the environment can choose to give you lower threads.
So we have to check how many threads we actually got and how many threads we asked.What if the user asks for 40000 threads or
any large number of threads ? You cant ask for more threads.There is a max level upto which the system can allot you the threads
Anything above that is not possible for the system to allot.So if we ask for above the optimal level,then the system allocates 
the maximum possible and not what the user asked for.So we can cross check using the master thread.By master thread you can 
get the actual no of threads allocated by the system.

if you notice  the results,by increasing threads,there is a small issue you will find.The issue is that by increassing the 
no of threads,the computational time will fall first to some extent and after that it almost becomes constant.Optimising this
futher is a bit difficult.There is also one more issue know as FALSE SHARING which we havent dealt with in our program.
This issue can be solved using synchronisation.

*/

/*program-4:- In the above problem,you probably used an array to create space for each thread to store its partial sum.If array
elements happen to store a cache line,this leads to false sharing.Modify your above program to avoid false sharing due to the 
sum array using synchronisation */


int main()
{
	int num_steps=100000;
	double step=1.0/num_steps;
	int num_threads=3;
	omp_set_num_threads(num_threads);
	double pi=0.0;
	
	#pragma omp parallel
	{
		double x;
		double sum=0.0;
		int i;
		for(i=0;i<num_steps;i=i+num_threads)
		{
		    x=(i+0.5)*step;
			sum+=4.0/(1+x*x);	
		}
		#pragma omp critical
		{
			pi+=sum*step;
		}
	}
	cout<<pi<<endl;
}

/*there is also one more issue we have to deal with.Sometimes its difficult to spilt the loops if we have loop dependent 
varaibles inside them. So we have to make sure that we dont have any loop dependent varibales inside our loops to make the 
splitting easy.Lets take one example below and understand this carefully.*/

/*lets say we want to compute the average value of an array of elements.*/
int MAX = 100000;
double ave = 0.0;
double A[MAX];
for(int i=0;i<MAX;i++)
{
   ave = ave +A[i]; //computing the sum of all the elements of A and storing in ave variable
}
ave = ave/MAX; //we will now divide by the number of elements in the array to get the average 

/*notice the above programe carefully.it has loop dependencies.The value of ave variable in loop i,is dependent on the value
of the ave variable in loop i-1. So we have loop dependencies and in this case its difficult to the spilt the work among
different threads. So we have to take care of these slight things and make the program as loop independent as possible.








