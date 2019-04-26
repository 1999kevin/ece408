#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "import.h"

void matrix_T(const float **a_matrix, float **a_t, int row, int col)
{
	int i, j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			a_t[j][i] = a_matrix[i][j];
		}
	}
}

void coo_T(const sparse_coo *a_coo, sparse_coo *a_t)
{
	for (int i = 0; i < a_coo->N; i++) {
		a_t->v[i] = a_coo->v[i];
		a_t->i[i] = a_coo->j[i];
		a_t->j[i] = a_coo->i[i];
	}
}

void matrix_m1(const float **a_matrix, const float **b_matrix, float **c_matrix,
int A_row, int A_col, int B_row, int B_col)
{
	if (A_row != B_col || A_col != B_row) {
		return;
	}
	
	int i, j, k;
	
	for (i = 0; i < A_row; i++)
	{
		for (j = 0; j < B_col; j++)
		{
			c_matrix[i][j] = 0;
			for (k = 0; k < A_col; k++) {
				c_matrix[i][j] += a_matrix[i][k] * b_matrix[k][j];
			}
		}
	}
}


// find (row,col) in A, if it's nonzero, return true and let e=(row,col)
int matrix_find(const sparse_coo *A, int row, int col, float *e)
{
	int p;
	for (p=0;p<A->N;p++){
		if (A->i[p]==row && A->j[p]==col){
			*e = A->v[p];
			return 0;
		}
	}
	
	return -1;
}

// coo-matrix*matrix
sparse_coo * matrix_m3(const sparse_coo *A, const sparse_coo *B, int L_1)
{
	// csr_B.T
	sparse_coo *C;
	int i,j,k,p;
	float m,t,s;
	int *i_tmp, *j_tmp;
	float *v_tmp;
	int elem_max;
	if (A->N > B->N) {
		elem_max = A->N;
	} else {
		elem_max = B->N;
	}
	i_tmp = (int *) malloc(sizeof(int)*elem_max);
	j_tmp = (int *) malloc(sizeof(int)*elem_max);
	v_tmp = (float *) malloc(sizeof(float)*elem_max);
	
	p = 0;
	for (i = 0; i < A->m; i++) {			// loop C-row
		for (j = 0; j < B->n; j++) {		// loop C-col
			s = 0;
			for (k = 0; k < A->n; k++) {	// loop A-col
				if (-1==matrix_find(A,i,k,&m)) {		// determine A[i,k] = 0?
					continue;
				}
				if (-1==matrix_find(B,k,j,&t)) {		// determine B[k,j] = 0?
					continue;
				}
				s = s + m*t;
			}
			
			if (s != 0) {		// C[i,j]!=0
				i_tmp[p] = i;
				j_tmp[p] = j;
				v_tmp[p] = s;
				p++;
			}
		}
	}

	C = sparse_coo_create(A->m, B->n, p);

	for (i = 0; i < p; i++) {
		C->v[i] = v_tmp[i] / L_1;
		C->i[i] = i_tmp[i];
		C->j[i] = j_tmp[i];
	}

	free(i_tmp);
	free(j_tmp);
	free(v_tmp);

	return C;
}

sparse_coo * toe_full_elem(const float *row0, const float **A, int n)
{
	sparse_coo *Out;
	int i, j, count;
	int *i_tmp, *j_tmp;
	float *v_tmp;
	i_tmp = (int *) malloc(sizeof(int)*n*n);
	j_tmp = (int *) malloc(sizeof(int)*n*n);
	v_tmp = (float *) malloc(sizeof(float)*n*n);
	count = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (row0[abs(j-i)] != 0 && A[i][j] != 0) {
				v_tmp[count] = row0[abs(j-i)] * A[i][j];
				i_tmp[count] = i;
				j_tmp[count] = j;
				count++;
			}
		}
	}

	Out = sparse_coo_create(n, n, count);

	for (i = 0; i < count; i++) {
		Out->v[i] = v_tmp[i];
		Out->i[i] = i_tmp[i];
		Out->j[i] = j_tmp[i];
	}

	free(i_tmp);
	free(j_tmp);
	free(v_tmp);
	return Out;
}

void matrix_trans(const sparse_coo *A, sparse_rcs *B) {
	int p;
	for (p = 0; p < B->N; p++) {
		B->v[p] = A->v[p];
		B->j[p] = A->j[p];
	}
	
	int t = 0; // counter for # nonzero data
	p = 0;
	while (p < B->m) {
		if (A->i[t] > p) {
			printf("%d  %d\n",p,A->i[t]);
			B->r[p] = t;
			p++;
		}
		else if (A->i[t] == p) {
			//printf("%d  %d\n",p,A->i[t]);
			B->r[p] = t;
			p++;
			t++;
		}
		else{
			t++;
		}
	}
	B->r[p] = B->N;
}
