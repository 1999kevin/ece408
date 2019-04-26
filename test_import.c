#include "import.h"

int main(int argc, char **argv){
  sparse_coo *H;
  printf("import H\n");
  H = sparse_coo_import("H.csv");
  printf("Finish importing H\n");
  
  full_r *e;
  printf("import e\n");
  e = full_r_import("e.csv");
  printf("Finish importing e\n");

  toeplitz *C;
  printf("import C\n");
  C = toeplitz_import("C.csv");
  printf("Finish importing C\n");
  return 0;
}
