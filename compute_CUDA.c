#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <string.h>
#include <assert.h>

#include "import.h"
#include "matrix.h"
#include "type.h"

void denseMatrixMultiply_wrapper(float **H, float **B, float **C, int numARows, int numAColumns, int numBRows, int numBColumns);
void denseMatrixTranspose_wrapper(float **H, float **A_T, int numARows, int numAColumns);
void cooToFullRank_wrapper(sparse_coo *A_coo, float **A_FR);
void sparseVector_wrapper(sparse_rcs *H, float **B, float **C, int FR_RC, int Div);

int main() {
	int i;
	sparse_rcs *H = NULL;
	printf("import H\n");
	FILE *fid;
	int m, n, N;

	assert("data/H.csv");

	fid = fopen("data/H.csv", "r");
	assert(fid);

	/* Read csv */
	char first_line[100];
	char *m_c, *n_c, *N_c;
	char *line;
	char *p1;
	int len;

	char *p_end;
	char *data_line;

	char tmp[25];

	fgets(first_line, 100, fid);
	line = (char *)first_line;

	for (int i = 0; i < 2; i++) {
		p1 = strchr(line, ',');
		len = strlen(line) - strlen(p1);
		switch (i) {
		case 0:
			N_c = (char *)malloc(len + 1);
			strncpy(N_c, line, len);
			N_c[len] = '\0';
			break;
		case 1:
			m_c = (char *)malloc(len + 1);
			strncpy(m_c, line, len);
			m_c[len] = '\0';
			break;
		}
		line = p1 + 1;
	}

	n_c = (char *)malloc(strlen(line));
	strncpy(n_c, line, strlen(line));

	N = (int)strtof(N_c, &p_end);
	m = (int)strtof(m_c, &p_end);
	n = (int)strtof(n_c, &p_end);

	H = sparse_rcs_create(m, n, N);

	/* Data vector */
	data_line = (char *)malloc(25 * N);
	fgets(data_line, 25 * N, fid);

	line = data_line;
	for (i = 0; i < N - 1; i++) {
		p1 = strchr(line, ',');
		len = strlen(line) - strlen(p1);
		strncpy(tmp, line, len);
		tmp[len] = '\0';
		line = p1 + 1;
		(H->v)[i] = strtof(tmp, &p_end);
	}
	(H->v)[i] = strtof(line, &p_end);

	/* Col vector */
	fgets(data_line, 25 * N, fid);

	line = data_line;
	for (i = 0; i < N - 1; i++) {
		p1 = strchr(line, ',');
		len = strlen(line) - strlen(p1);
		strncpy(tmp, line, len);
		tmp[len] = '\0';
		line = p1 + 1;
		(H->j)[i] = (int)strtof(tmp, &p_end);
	}
	(H->j)[i] = (int)strtof(line, &p_end);

	/* PTR */
	fgets(data_line, 25 * (m + 1), fid);

	line = data_line;
	for (i = 0; i < m; i++) {
		p1 = strchr(line, ',');
		len = strlen(line) - strlen(p1);
		strncpy(tmp, line, len);
		tmp[len] = '\0';
		line = p1 + 1;
		(H->r)[i] = (int)strtof(tmp, &p_end);
	}
	(H->r)[i] = (int)strtof(line, &p_end);

	free(N_c);
	free(m_c);
	free(n_c);
	free(data_line);

	fclose(fid);

	//H = sparse_rcs_import("data/H.csv");
	printf("Finish importing H\n");
	full_r *e;
	printf("import e\n");
	e = full_r_import("data/e.csv");
	printf("Finish importing e\n");
	toeplitz *C;
	printf("import C\n");
	C = toeplitz_import("data/C.csv");
	printf("Finish importing C\n");


	/* Compute e @ e.T */
	float **e_2;
	float **e_T;
	e_T = malloc(sizeof(float*)*e->n);
	for (i = 0; i < e->n; i++) {
		e_T[i] = malloc(sizeof(float)*e->m);
	}
	e_2 = malloc(sizeof(float*)*e->m);
	for (i = 0; i < e->m; i++) {
		e_2[i] = malloc(sizeof(float)*e->m);
	}
	denseMatrixTranspose_wrapper(e->v, e_T, e->m, e->n);
	//matrix_T((const float **) e->v, e_T, e->m, e->n);
	printf("Transpose e\n");
	denseMatrixMultiply_wrapper(e->v, e_T, e_2, e->m, e->n, e->n, e->m);
	//matrix_m1((const float **) e->v, (const float **) e_T, e_2, e->m, e->n, e->n, e->m);
	printf("Compute e^2\n");

	/* Compute np.multiply(C, e @ e.T) */
	sparse_coo *C_e_2;
	printf("Compute np.multiply(C, e @ e.T)\n");
	C_e_2 = toe_full_elem((const float *) C->row0, (const float **) e_2, C->n);
	printf("Num nonzero elements of C * e^2: %d\n", C_e_2->N);

	/* C_e_2 -> full rank */
	float **C_e_2_FR;
	C_e_2_FR = malloc(sizeof(float*)*C_e_2->n);
	for (i = 0; i < e->m; i++) {
		C_e_2_FR[i] = malloc(sizeof(float)*C_e_2->m);
	}
	cooToFullRank_wrapper(C_e_2, C_e_2_FR);
	printf("C_e_2 Full Rank\n");

	/* Compute np.multiply(C, e @ e.T) @ H.T / (L-1) */
	/* (H @ C_e_2.T).T */
	float **P_HT_FR;
	P_HT_FR = malloc(sizeof(float*)*C_e_2->n);
	for (i = 0; i < C_e_2->m; i++) {
		P_HT_FR[i] = malloc(sizeof(float)*(H->m));
	}
	sparseVector_wrapper(H, C_e_2_FR, P_HT_FR, C_e_2->n, e->n - 1);

	/* Export P_HT to .csv file in CSR format */
	printf("Write P_HT to .csv file in CSR format\n");

	//sparse_rcs *P_HT;
	int j, p, k;
	p = 0;
	//calculate N
	for (i = 0; i<C_e_2->n; i++) {
		for (j = 0; j<H->m; j++) {
			if (P_HT_FR[i][j] != 0) {
				//printf("%f i:%d j:%d\n",A[i][j],i,j);
				p++;
			}
		}
	}
	//printf("wotaicaile1\n");
	//create P_HT
	sparse_rcs* P_HT = sparse_rcs_create(C_e_2->n, H->m, p);
	//determine v,j,r
	p = 0;
	//initialize
	P_HT->r[0] = 0;
	for (i = 0; i<C_e_2->n; i++) {
		for (j = 0; j<H->m; j++) {
			if (P_HT_FR[i][j] != 0 && p<P_HT->N) {
				P_HT->v[p] = P_HT_FR[i][j];
				P_HT->j[p] = j;

				p++;
			}
		}

		P_HT->r[i + 1] = p;
	}
	fid = fopen("data/P_HT.csv", "w");
	fprintf(fid,"%d,",P_HT->N);
	fprintf(fid,"%d,",P_HT->m);
	fprintf(fid,"%d\n",P_HT->n);
	for (i = 0; i < P_HT->N-1; i++) {
		fprintf(fid,"%f,",(P_HT->v)[i]);
	}
	fprintf(fid,"%f\n",(P_HT->v)[i]);
	for (i = 0; i < P_HT->N-1; i++) {
		fprintf(fid,"%d,",(P_HT->j)[i]);
	}
	fprintf(fid,"%d\n",(P_HT->j)[i]);
	for (i = 0; i < P_HT->m; i++) {
		fprintf(fid,"%d,",(P_HT->r)[i]);
	}
	fprintf(fid,"%d\n",(P_HT->r)[i]);
	fclose(fid);

	/* Free all structures */
	sparse_rcs_free(H);
	sparse_coo_free(C_e_2);

	//sparse_rcs_free(P_HT);

	full_r_free(e);

	toeplitz_free(C);

	for (i = 0; i < e->n; i++) {
		free(e_T[i]);
	}
	free(e_T);
	for (i = 0; i < e->m; i++) {
		free(e_2[i]);
	}
	free(e_2);
	for (i = 0; i < C->n; i++) {
		free(C_e_2_FR[i]);
	}
	free(C_e_2_FR);
	for (i = 0; i < C->n; i++) {
		free(P_HT_FR[i]);
	}
	free(P_HT_FR);
	printf("Free all the memory\n");

	return 0;
}
