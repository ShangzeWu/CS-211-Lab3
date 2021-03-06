/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   int   local_first;
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   char  *local_prime_marked;
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;
   unsigned long int  local_prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */
   unsigned long int  local_prime_size;
   char  *marked2; 

   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);
   
   

   unsigned long int oddn = n - n / 2 - 1;  //排除1和偶数之后的总数量

   unsigned long int size2 = (int) sqrt((double) n) + 1; //根号n的大小
                     size2 = size2 - size2 / 2 - 1;      //去除偶数部分和1

   unsigned long int low_value_idx = id * oddn / p; //同一个processor中最小的下标
   unsigned long int high_value_idx = -1 + (id + 1) * oddn / p; //最大的下标
   size = high_value_idx - low_value_idx + 1; //array的大小
   low_value = 2 * low_value_idx + 3;  //array中最小的值
   high_value = 2 * high_value_idx + 3; //array中最大的值

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (oddn - 1) / p;

   if (proc0_size < (int) sqrt((double) oddn)) {
      if (!id) printf("Too many processes\n");
      MPI_Finalize();
      exit(1);
   }


   marked = (char *) malloc(size);
   marked2 = (char *) malloc(size2); 

   if (marked == NULL) {
      printf("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }

   for (i = 0; i < size; i++) marked[i] = 0;  //marked 全部赋值成0
   for (i = 0; i < size2; i++) marked2[i] = 0; 

   //以下do-while循环用于在各个processor中计算根号n范围内的质数，如果是非质数，标注为1.
   index = 0;
   prime = 3;
   int tempn = (int) sqrt((double) n);
   do {
      for (i = (prime * 3 - 3) / 2; i < size2; i += prime) marked2[i] = 1;
         while (marked2[++index]!=0);
         prime = 2 * index + 3;
   } while (prime * prime <= tempn);


   unsigned long int blockSize = 100000;
   unsigned long int itr = size / blockSize;
   unsigned long int blockStart = 0;
   if (size % blockSize) itr++; //除不尽的情况下， itr加1

   while(itr--)
   {
      low_value = 2 * (low_value_idx + blockStart) + 3;
      high_value = 2 * (low_value_idx + blockStart + blockSize) + 3;
      if (itr == 1)
         high_value = 2 * high_value_idx + 3;
      index = 0;
      prime = 3;
      do {
         if (prime * prime > low_value)   //筛选processor， 0号满足，1号及以后不满足
            first = (prime * prime - low_value) / 2;
         //odd number minuses odd number = even number, divide by 2 to find index
         else {   //1号及以后不满足
            if ((low_value % prime)==0) first = 0;
            else 
            {
               first = (low_value / prime + 1) * prime;
               first = ((first - low_value) % 2) == 0 ? first : first + prime; //确保是奇数
               first = (first - low_value) / 2;  //换算成对应的array中的index
            }
         }
         for (i = blockStart + first; i < blockStart + blockSize && i < size; i += prime) marked[i] = 1;
      
         while (marked2[++index]!=0); //循环执行到当前maeked[index]==0， 保留index的值
         prime = 2 * index + 3;   //将index转换成对应的prime。

      } while (prime * prime <= high_value); //停止寻找下一个prime。
      
      blockStart += blockSize;
   }


   count = 0;
   for (i = 0; i < size; i++)
      if (marked[i]==0) count++; //统计当前processor下的prime数量
   if (p > 1)
      MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                  0, MPI_COMM_WORLD);

   global_count++;
   //算上‘2’ 


   elapsed_time += MPI_Wtime();

   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}

