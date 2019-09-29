#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


int main(int argc, char* argv[]) {

  int x = 100, i, n = 100;
  //double a[(int)sqrt(x)];
  double a[n];
  double b[n];
  
  printf("%d\n", (int)sqrt(x));
  int dotp = 0;
# pragma omp parallel for  
  /*for(i = 0; i < n; i++)
    dotp += a[i] * b[i];*/
  for (i = 0; i < n; i++) {
    a[i] = i*3;
    //if(a[i] < b[i]) break;
  }
  /*for(i = 0; i < (int)sqrt(x); i++)
{
  a[i] = 2.3 * i;
  printf("thread: %d, a[i]: %f\n", omp_get_thread_num(), a[i]);
  if(i < 10) b[i] = a[i];
}
  */

 for (i = 0; i < (int)sqrt(x); i++) {
   printf("%f\n", a[i]);
 }
 printf("\n");
}
