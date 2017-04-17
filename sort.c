#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

#include "mpi.h"

#define SIZE 10



static int intcompare(const void *i, const void *j)
{
  if ((*(int *)i) > (*(int *)j))
    return (1);
  if ((*(int *)i) < (*(int *)j))
    return (-1);
  return (0);
}



 main (int argc, char *argv[])
{
  

  int 	     Numprocs,MyRank, Root = 0;
  int 	     i,j,k, NoofElements, NoofElements_Bloc,
				  NoElementsToSort;
  int 	     count, temp, argument1;
  int 	     *Input, *InputData;
  int 	     *Splitter, *AllSplitter;
  int 	     *Buckets, *BucketBuffer, *LocalBucket;
  int 	     *OutputBuffer, *Output;
  FILE 	     *InputFile, *fp;
  MPI_Status  status; 
  
  
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

  if(argc != 2) {
      if(MyRank ==0) printf(" Usage : run size\n");
	 MPI_Finalize();
	 exit(0);
	}
  
  MPI_Barrier(MPI_COMM_WORLD);
  timestamp_type time1, time2;
  get_timestamp(&time1);

  
  
  if (MyRank == Root){
    argument1 =  atoi(argv[1]);
    NoofElements = argument1*Numprocs; 
    //NoofElements = atoi(argv[1]);
    Input = (int *) malloc (NoofElements*sizeof(int));
	 if(Input == NULL) {
		printf("Error : Can not allocate memory \n");
    }

  
    srand48((unsigned int)NoofElements);
	 for(i=0; i< NoofElements; i++) {
       Input[i] = rand();
       
    }
  }
  
  MPI_Bcast (&NoofElements, 1, MPI_INT, 0, MPI_COMM_WORLD);
 

  NoofElements_Bloc = NoofElements / Numprocs;
  InputData = (int *) malloc (NoofElements_Bloc * sizeof (int));

  MPI_Scatter(Input, NoofElements_Bloc, MPI_INT, InputData, 
				  NoofElements_Bloc, MPI_INT, Root, MPI_COMM_WORLD);

  
  qsort ((char *) InputData, NoofElements_Bloc, sizeof(int), intcompare);

  
  Splitter = (int *) malloc (sizeof (int) * (Numprocs-1));
  for (i=0; i< (Numprocs-1); i++){
        Splitter[i] = InputData[NoofElements/(Numprocs*Numprocs) * (i+1)];
  } 

  
  AllSplitter = (int *) malloc (sizeof (int) * Numprocs * (Numprocs-1));
  MPI_Gather (Splitter, Numprocs-1, MPI_INT, AllSplitter, Numprocs-1, 
				  MPI_INT, Root, MPI_COMM_WORLD);

  
  if (MyRank == Root){
    qsort ((char *) AllSplitter, Numprocs*(Numprocs-1), sizeof(int), intcompare);

    for (i=0; i<Numprocs-1; i++)
      Splitter[i] = AllSplitter[(Numprocs-1)*(i+1)];
  }
  
  
  MPI_Bcast (Splitter, Numprocs-1, MPI_INT, 0, MPI_COMM_WORLD);

  
  Buckets = (int *) malloc (sizeof (int) * (NoofElements + Numprocs));
  
  j = 0;
  k = 1;

  for (i=0; i<NoofElements_Bloc; i++){
    if(j < (Numprocs-1)){
       if (InputData[i] < Splitter[j]) 
			 Buckets[((NoofElements_Bloc + 1) * j) + k++] = InputData[i]; 
       else{
	       Buckets[(NoofElements_Bloc + 1) * j] = k-1;
		    k=1;
			 j++;
		    i--;
       }
    }
    else 
       Buckets[((NoofElements_Bloc + 1) * j) + k++] = InputData[i];
  }
  Buckets[(NoofElements_Bloc + 1) * j] = k - 1;
      
  

  BucketBuffer = (int *) malloc (sizeof (int) * (NoofElements + Numprocs));

  MPI_Alltoall (Buckets, NoofElements_Bloc + 1, MPI_INT, BucketBuffer, 
					 NoofElements_Bloc + 1, MPI_INT, MPI_COMM_WORLD);

  
  LocalBucket = (int *) malloc (sizeof (int) * 2 * NoofElements / Numprocs);

  count = 1;

  for (j=0; j<Numprocs; j++) {
  k = 1;
    for (i=0; i<BucketBuffer[(NoofElements/Numprocs + 1) * j]; i++) 
      LocalBucket[count++] = BucketBuffer[(NoofElements/Numprocs + 1) * j + k++];
  }
  LocalBucket[0] = count-1;
    
 

  NoElementsToSort = LocalBucket[0];
  qsort ((char *) &LocalBucket[1], NoElementsToSort, sizeof(int), intcompare); 

  
  if(MyRank == Root) {
  		OutputBuffer = (int *) malloc (sizeof(int) * 2 * NoofElements);
  		Output = (int *) malloc (sizeof (int) * NoofElements);
  }

  MPI_Gather (LocalBucket, 2*NoofElements_Bloc, MPI_INT, OutputBuffer, 
				  2*NoofElements_Bloc, MPI_INT, Root, MPI_COMM_WORLD);

  
	if (MyRank == Root){
		count = 0;
		for(j=0; j<Numprocs; j++){
          k = 1;
      	 for(i=0; i<OutputBuffer[(2 * NoofElements/Numprocs) * j]; i++) 
				 Output[count++] = OutputBuffer[(2*NoofElements/Numprocs) * j + k++];
    	}

      
    	if ((fp = fopen("sort.out", "w")) == NULL){
         	printf("Can't Open Output File \n");
      		exit(0);
    	}
		 
    	fprintf (fp, "Number of Elements to be sorted : %d \n", NoofElements);
    	
    	fprintf (fp, "The sorted sequence is : \n");
	
    	for (i=0; i<NoofElements; i++){
	      	fprintf(fp, "%d\n", Output[i]);
	      	//printf( "%d   ", Output[i]);
	}
	
	fclose(fp);
    	free(Input);
  	free(OutputBuffer);
  	free(Output);
   }

  	free(InputData);
  	free(Splitter);
  	free(AllSplitter);
  	free(Buckets);
  	free(BucketBuffer);
  	free(LocalBucket);

  /**timming**/
  MPI_Barrier(MPI_COMM_WORLD);
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  if (0 == MyRank) {
    printf("Time elapsed is %f seconds.\n", elapsed);
    fprintf (fp, "time Elapsed %f \n", elapsed);
  }

   /**** Finalize ****/
   MPI_Finalize();
}
