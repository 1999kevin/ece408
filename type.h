#ifndef TYPE_H
#define TYPE_H

typedef struct {
  int N;	// non-zero num
  int m, n;	// row, col
  float *v;	// data vec
  int *j;	// col vec
  int *i;	// row vec
} sparse_coo;

typedef struct {
	int N;
	int m, n;
	float *v;
	int *j;
	int *r;
} sparse_rcs;

typedef struct {
	int m, n;
	float **v;
} full_r;

typedef struct {
	int N, n;
	float *row0;
} toeplitz;

#endif
