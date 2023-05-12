#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mpi.h>

void memtrack(char *message, int *myrank)

{

        struct rusage myrusage;



         //MPI_Barrier(MPI_COMM_WORLD); // optional, depending on your use case

                getrusage(RUSAGE_SELF, &myrusage);

                         printf("%d: %s: maxrss=%.1fMB\n",

                                        *myrank, message, myrusage.ru_maxrss/1024.0);

                                         }



                                          void  memtrack_(int *myrank)

                                         {

                                                 memtrack("", myrank);

                                                 }
