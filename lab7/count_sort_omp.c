/*

1. private(count, i, j), shared(a, n)
2. No, each current iteration only relies on the values of a at i and j, and count and each count will be different.
3. 
for (i = 0; i < n; i++) {                                                                
   a[i] = temp[i];                                                                        
}

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

int cmpfunc (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    return -1;
  }
  int num_threads;
  num_threads = strtol(argv[1], NULL, 10);
  int i, j, count, k;
  //int* temp = (int*) malloc(n*sizeof(int));  
  int n = 10000;
  int a[n];
  int p[] = {1, 2, 4, 6, 8, 10, 12, 14, 16};
  for (i = 0; i < n; i++) {
    a[i] = n - i;
  }
  int* temp = (int*) malloc(n*sizeof(int));
  double totalTime;
  //for (k = 0; k < 9; k++) {
  totalTime = omp_get_wtime();
# pragma omp parallel for num_threads(num_threads) private(count, j, i) shared(a, n)
for(i=0; i<n; i++)
{
  //printf("%d\n", omp_get_thread_num());
  count = 0;
  for(j=0; j<n; j++)
  {
    if(a[j]<a[i]) count++;
    else if(a[j]==a[i] && j<i) count++;
  }
  # pragma omp critical
  temp[count] = a[i];
}
for (i = 0; i < n; i++) {                                                               
  a[i] = temp[i];                                                                        
}
printf("parallel  time: %f\n", k, omp_get_wtime() - totalTime);
//totalTime = omp_get_wtime();
//qsort(a, n, sizeof(int), cmpfunc);
//printf("qsort  time: %f\n", k, omp_get_wtime() - totalTime);
 /*for (i = 0; i < n; i++) {
   a[i] = temp[i];
 }
 printf("parallel  time: &f\n", k, omp_get_wtime() - totalTime);
 */
 // serial
 /* totalTime = omp_get_wtime();
for(i=0; i<n; i++)
{
  //printf("%d\n", omp_get_thread_num());                                             
  count = 0;
  for(j=0; j<n; j++)
  {
    if(a[j]<a[i]) count++;
    else if(a[j]==a[i] && j<i) count++;
  }
  temp[count] = a[i];
}
 printf("serial time: %f\n", omp_get_wtime() - totalTime); 
 */
/*for (i = 0; i < n; i++) {
   a[i] = temp[i];
 }
 printf("serial time: &f\n", omp_get_wtime() - totalTime);
*/
free(temp);


/*for (i = 0; i < n; i++) {
   printf("%d: %d\n", i, a[i]);
 }*/
}
