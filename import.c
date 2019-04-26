#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "import.h"
//#include "blas.h"
//#include "util.h"

sparse_rcs *sparse_rcs_create(int m, int n, int N) {
  sparse_rcs *A;

  assert(m >= 0);
  assert(n >= 0);
  assert(N >= 0);

  A = malloc(sizeof(sparse_rcs));
  assert(A);

  A->m = m;
  A->n = n;
  A->N = N;

  if (N == 0) {
    A->v = NULL;
    A->j = NULL;
    A->r = NULL;
  }
  else {
    A->v = malloc(sizeof(float) * N);
    assert(A->v);

    A->j = malloc(sizeof(float) * N);
    assert(A->j);

    A->r = malloc(sizeof(float) * (m+1));
    assert(A->r);
  }

  return A;
}

void sparse_rcs_free(sparse_rcs *A) {
  free(A->v);
  free(A->j);
  free(A->r);
  free(A);
}

sparse_coo *sparse_coo_create(int m, int n, int N) {
  sparse_coo *A;

  assert(m >= 0);
  assert(n >= 0);
  assert(N >= 0);

  A = malloc(sizeof(sparse_coo));
  assert(A);

  A->m = m;
  A->n = n;
  A->N = N;

  if (N == 0) {
    A->v = NULL;
    A->j = NULL;
    A->i = NULL;
  }
  else {
    A->v = malloc(sizeof(float) * N);
    assert(A->v);

    A->j = malloc(sizeof(float) * N);
    assert(A->j);

    A->i = malloc(sizeof(float) * N);
    assert(A->i);
  }

  return A;
}

void sparse_coo_free(sparse_coo *A)
{
  free(A->v);
  free(A->i);
  free(A->j);
  free(A);
}

sparse_coo *sparse_coo_import(const char *filename){
  FILE *fid;
  sparse_coo *A = NULL;
  int m, n, N;

  assert(filename);

  fid = fopen(filename, "r");
  assert(fid);

  /* Read csv */
  char first_line[100];
  char *m_c, *n_c, *N_c;
  char *line;
  char *p1;
  int len;

  char *p_end;
  char *data_line;

  char tmp[20];
  int i;

  fgets(first_line, 100, fid);
  line = (char *)first_line;

  for (int i = 0; i < 2; i++) {
    p1 = strchr(line,',');
    len = strlen(line) - strlen(p1);
    switch (i) {
      case 0:
        N_c = (char *) malloc(len+1);
        strncpy(N_c, line, len);
        N_c[len] = '\0';
        break;
      case 1:
        m_c = (char *) malloc(len+1);
        strncpy(m_c, line, len);
        m_c[len] = '\0';
        break;
    }
    line = p1 + 1;
  }

  n_c = (char *) malloc(strlen(line));
  strncpy(n_c, line, strlen(line));

  N = (int) strtof(N_c, &p_end);
  //printf("N: %d\n", N);
  m = (int) strtof(m_c, &p_end);
  //printf("m: %d\n", m);
  n = (int) strtof(n_c, &p_end);
  //printf("n: %d\n", n);

  A = sparse_coo_create(m, n, N);

  /* Data vector */
  data_line = (char *) malloc(21*N);
  fgets(data_line, 21*N, fid);

  line = data_line;
  for (i = 0; i < N-1; i++) {
    p1 = strchr(line,',');
    len = strlen(line) - strlen(p1);
    strncpy(tmp, line, len);
    tmp[len] = '\0';
    line = p1 + 1;
    (A->v)[i] = strtof(tmp, &p_end);
  }
  (A->v)[i] = strtof(line, &p_end);
  /*printf("Data vector:\n" );
  for (i = 0; i < N; i++)
    printf("%f ", (A->v)[i]);
  printf("\n");*/

  /* Col vector */
  fgets(data_line, 21*N, fid);

  line = data_line;
  for (i = 0; i < N-1; i++) {
    p1 = strchr(line,',');
    len = strlen(line) - strlen(p1);
    strncpy(tmp, line, len);
    tmp[len] = '\0';
    line = p1 + 1;
    (A->j)[i] = (int) strtof(tmp, &p_end);
  }
  (A->j)[i] = (int) strtof(line, &p_end);
  /*printf("Col vector:\n" );
  for (i = 0; i < N; i++)
    printf("%d ", (A->j)[i]);
  printf("\n");*/

  /* Row vector */
  fgets(data_line, 21*N, fid);

  line = data_line;
  for (i = 0; i < N-1; i++) {
    p1 = strchr(line,',');
    len = strlen(line) - strlen(p1);
    strncpy(tmp, line, len);
    tmp[len] = '\0';
    line = p1 + 1;
    (A->i)[i] = (int) strtof(tmp, &p_end);
  }
  (A->i)[i] = (int) strtof(line, &p_end);
  /*printf("Row vector:\n" );
  for (i = 0; i < N; i++)
    printf("%d ", (A->i)[i]);
  printf("\n");*/

  free(N_c);
  free(m_c);
  free(n_c);
  free(data_line);

  fclose(fid);

  return A;
}

full_r *full_r_create(int m, int n) {
  assert(m >= 0);
  assert(n >= 0);

  full_r *A;
  A = malloc(sizeof(full_r));
  assert(A);

  A->m = m;
  A->n = n;
  if (m > 0) {
    A->v = malloc(sizeof(float *)*m);
    assert(A->v);
  }

  for (int i = 0; i < m; i++) {
    (A->v)[i] = malloc(sizeof(float)*n);
  }

  return A;
}

void full_r_free(full_r *A)
{
  for (int i = 0; i < A->m; i++) {
    free(A->v[i]);
  }
  free(A->v);
  free(A);
}

full_r *full_r_import(const char *filename) {
  FILE *fid;
  full_r *A = NULL;
  int m, n;

  assert(filename);

  fid = fopen(filename, "r");
  assert(fid);

  char first_line[100];
  char *m_c;
  char *line;
  char *p1;
  int len;

  char *p_end;
  char *data_line;

  char tmp[20];
  int row, i;

  fgets(first_line, 100, fid);
  line = (char *)first_line;

  p1 = strchr(line,',');
  len = strlen(line) - strlen(p1);
  m_c = (char *) malloc(len+1);
  strncpy(m_c, line, len);
  m_c[len] = '\0';
  line = p1 + 1;

  m = (int) strtof(m_c, &p_end);
  //printf("m: %d\n", m);
  n = (int) strtof(line, &p_end);

  A = full_r_create(m, n);

  data_line = (char *) malloc(21*n);
  
  //printf("Data Matrix:\n");
  for (row = 0; row < m; row++) {
    fgets(data_line, 21 * n, fid);
    line = data_line;
    for (i = 0; i < n-1; i++) {
      p1 = strchr(line,',');
      len = strlen(line) - strlen(p1);
      strncpy(tmp, line, len);
      tmp[len] = '\0';
      line = p1 + 1;
      (A->v)[row][i] = strtof(tmp, &p_end);
    }
    (A->v)[row][i] = strtof(line, &p_end);
    /*for (i = 0; i < n; i++)
      printf("%f ", (A->v)[row][i]);
    printf("\n");*/
  }

  free(m_c);
  free(data_line);

  fclose(fid);
  return A;
}

toeplitz *toeplitz_create(int N, int n) {
  assert(N >= 0);
  assert(n >= 0);

  toeplitz *A;
  A = malloc(sizeof(toeplitz));
  assert(A);

  A->N = N;
  A->n = n;
  if (n > 0) {
    A->row0 = malloc(sizeof(float)*n);
    assert(A->row0);
  }

  return A;
}

void toeplitz_free(toeplitz *A)
{
  free(A->row0);
  free(A);
}

toeplitz *toeplitz_import(const char *filename) {
  FILE *fid;
  toeplitz *A = NULL;
  int n, N;

  assert(filename);

  fid = fopen(filename, "r");
  assert(fid);

  char first_line[100];
  char *n_c, *N_c;
  char *line;
  char *p1;
  int len;

  char *p_end;
  char *data_line;

  char tmp[20];
  int i;

  fgets(first_line, 100, fid);
  line = (char *)first_line;

  for (int i = 0; i < 3; i++) {
    p1 = strchr(line,',');
    len = strlen(line) - strlen(p1);
    if (i == 2) {
      N_c = (char *) malloc(len+1);
      strncpy(N_c, line, len);
      N_c[len] = '\0';
    }
    line = p1 + 1;
  }
  n_c = (char *) malloc(strlen(line));
  strncpy(n_c, line, strlen(line));

  N = (int) strtof(N_c, &p_end);
  //printf("N: %d\n", N);
  n = (int) strtof(n_c, &p_end);
  //printf("n: %d\n", n);

  A = toeplitz_create(N, n);
  /* Data vector */
  data_line = (char *) malloc(21*n);
  fgets(data_line, 21*n, fid);

  line = data_line;
  for (i = 0; i < n-1; i++) {
    p1 = strchr(line,',');
    len = strlen(line) - strlen(p1);
    strncpy(tmp, line, len);
    tmp[len] = '\0';
    line = p1 + 1;
    (A->row0)[i] = strtof(tmp, &p_end);
  }
  (A->row0)[i] = strtof(line, &p_end);

  free(N_c);
  free(n_c);
  free(data_line);

  fclose(fid);
  return A;

}
