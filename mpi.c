/* MPI-parallel Jacobi smoothing to solve -u''=f
 * Global vector has N unknowns, each processor works with its
 * part, which has lN = N/p unknowns.
 * Author: Georg Stadler
 * Co-Author: Chirag Singhal
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

#include "mpi.h"




int main(int argc, char * argv[])
{
  int mpirank, p,j,i ; 
  MPI_Status status, status1;
  int sqrtp, N, lN, iter, max_iters;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);

  sscanf(argv[1], "%d", &N);
  sscanf(argv[2], "%d", &max_iters);

  sqrtp = sqrt(p);
  lN = N / sqrtp;
  if (mpirank == 0 && (N % sqrtp != 0) ) {
    printf("Exiting. N must be a multiple of p\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  timestamp_type time1, time2;
  get_timestamp(&time1);

  
  double **lu = (double **)malloc((lN+2) * sizeof(double *));
  for (i=0; i < (lN+2); i++)
	lu[i] = (double *)malloc((lN+2) * sizeof(double));
  
  
  double **lunew = (double **)malloc((lN+2) * sizeof(double *));
  for (i=0; i < (lN+2); i++)
        lunew[i] = (double *)malloc((lN+2) * sizeof(double));


  
  double h = 1.0 / (N + 1);
  double hsq = h * h;
  double invhsq = 1./hsq;
  
  for (iter = 0; iter < max_iters; iter++) {

    for (i = 1; i <= lN; i++){
        for(j = 1; j <= lN; j++){
          lunew[i][j] =  (hsq + lu[i-1][j] + lu[i][j-1] + lu[i+1][j] + lu[i][j+1]) /4;
        }  
    }

        if(mpirank < p - sqrtp){
              for(i = 1; i<= lN; i++){
                      MPI_Send(&(lunew[i][lN]), 1, MPI_DOUBLE, mpirank+sqrtp, i, MPI_COMM_WORLD);
                              MPI_Recv(&(lunew[i][lN+1]), 1, MPI_DOUBLE, mpirank+sqrtp, i, MPI_COMM_WORLD, &status1);
                                    }
                                        }

   if(mpirank % sqrtp != 0){
      for(i = 1; i<= lN; i++){
        	MPI_Send(&(lunew[1][i]), 1, MPI_DOUBLE, mpirank-1, i, MPI_COMM_WORLD);
        	MPI_Recv(&(lunew[0][i]), 1, MPI_DOUBLE, mpirank-1, i, MPI_COMM_WORLD, &status1);
      }
    }


    if(mpirank % sqrtp != sqrtp -1){
      for(i = 1; i<= lN; i++){
        	MPI_Send(&(lunew[lN][i]), 1, MPI_DOUBLE, mpirank+1, i , MPI_COMM_WORLD);
        MPI_Recv(&(lunew[lN+1][i]), 1, MPI_DOUBLE, mpirank+1, i, MPI_COMM_WORLD, &status1); 
      }
    }
   
 

    if(mpirank >=  sqrtp){
      for(i = 1; i<= lN; i++){
        MPI_Send(&(lunew[i][1]), 1, MPI_DOUBLE, mpirank-sqrtp, i, MPI_COMM_WORLD);
       		 MPI_Recv(&(lunew[i][0]), 1, MPI_DOUBLE, mpirank-sqrtp, i, MPI_COMM_WORLD, &status1);
      }
    }
    

    
    for(i = 0; i < lN+2; i++) {
      for(j = 0; j < lN+2; j++){
        lu[i][j] = lunew[i][j];  
      }
    }
  }
  
  free(lu);
  free(lunew);

  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  if (0 == mpirank) {
    printf("Time elapsed is %f seconds.\n", elapsed);
  }
  MPI_Finalize();
  return 0;
}


double compute_residual(double *lu, int lN, double invhsq)
{
  int i;
  double tmp, gres = 0.0, lres = 0.0;

  for (i = 1; i <= lN; i++){
    tmp = ((2.0*lu[i] - lu[i-1] - lu[i+1]) * invhsq - 1);
    lres += tmp * tmp;
  }
  
  MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sqrt(gres);
}
