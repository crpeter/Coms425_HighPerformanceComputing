# include <stdio.h>

void print_option(int op, int n) {
  // printf("option: %d", op);
  int i, count;
  count = 0;
  if (op == 0) {
    for(i=0; i < n; ++i) {
      for(count = 1; count <= i+1; ++count) {
	printf("%d ", count);
      }
      printf("\n");
    }
  }
  else if (op == 1) {
    for(i=0; i < n; ++i) {
      for(count = 1; count <= n-i; ++count){
	printf("%d ", count);
      }
      printf("\n");
    }
  }
  else if (op == 2) {
    int helpCount = 1;
    for(i=0; i < n; ++i){
      helpCount = 1;
      for(count = 1; count <= n; ++count) {
	if (count < n-i ) {
	  printf("  ");
	}
	else {
	  printf("%d ", helpCount);
	  helpCount += 1;
	}
      }
      printf("\n");
    }
  }
  else if (op == 3) {
    for (i=0; i < n; ++i) {
      int helpCount = 1;
      for(count = 1; count <= n; ++count) {
	if (count <= i) {
	  printf("  ");
	} else {
	  printf("%d ", helpCount);
	  ++helpCount;
	}
      }
      printf("\n");
    }
  }
}

int main() {
  int c, n;
  printf("Enter Character:");
  scanf("%c", &c);
  printf("Enter Number:");
  scanf("%d", &n);
  //printf("\nSelected %d, %c", n, c);
  
  switch (c) {
  case 'A':
    print_option(0, n);
    break;
  case 'B':
    print_option(1, n);
    break;
  case 'C':
    print_option(2, n);
    break;
  case 'D':
    print_option(3, n);
    break;
  }
}
