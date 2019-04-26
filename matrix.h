#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include "type.h"

void matrix_T(const float **a_matrix, float **a_t, int row, int col);
void coo_T(const sparse_coo *a_coo, sparse_coo *a_t);
void matrix_m1(const float **a_matrix, const float **b_matrix, float **c_matrix, int A_row, int A_col, int B_row, int B_col);
int matrix_find(const sparse_coo *A, int row, int col, float *e);
sparse_coo * matrix_m3(const sparse_coo *A, const sparse_coo *B, int L_1);
sparse_coo * toe_full_elem(const float *row0, const float **A, int n);
void matrix_trans(const sparse_coo *A, sparse_rcs *B);
#endif
