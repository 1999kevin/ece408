#include <stdlib.h>
#include <stdio.h>

#include "import.h"
#include "matrix.h"

int main() {
	int i;
	sparse_coo *H;
	printf("import H\n");
	H = sparse_coo_import("data/H.csv");
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
	matrix_T((const float **) e->v, e_T, e->m, e->n);
	printf("Transpose e\n");
	/*for (i = 0; i < e->n; i++) {
		for (j = 0; j < e->m; j++) {
			printf("%f ", e_T[i][j]);
		}
		printf("\n");
		break;
	}*/
	matrix_m1((const float **) e->v, (const float **) e_T, e_2, e->m, e->n, e->n, e->m);
	printf("Compute e^2\n");
	/*for (i = 0; i < e->m; i++) {
		for (j = 0; j < e->m; j++) {
			printf("%f ", e_2[i][j]);
		}
		printf("\n");
		break;
	}*/

	/* Compute np.multiply(C, e @ e.T) */
	sparse_coo *C_e_2;
	printf("Compute np.multiply(C, e @ e.T)\n");
	C_e_2 = toe_full_elem((const float *) C->row0, (const float **) e_2, C->n);
	printf("Num nonzero elements of C * e^2: %d\n", C_e_2->N);
	/*for (i = 0; i < C_e_2->N; i++) {
		printf("%f ", C_e_2->v[i]);
	}
	printf("\n");*/

	/* Compute np.multiply(C, e @ e.T) @ H.T / (L-1) */
	sparse_coo *P_HT_coo;
	sparse_coo *H_T;
	H_T = sparse_coo_create(H->n, H->m, H->N);
	coo_T((const sparse_coo *) H, H_T);
	printf("Transpose H\n");
	P_HT_coo = matrix_m3((const sparse_coo *) C_e_2, (const sparse_coo *) H_T, e->n-1);
	printf("Compute np.multiply(C, e @ e.T) @ H.T\n");
	printf("Num nonzero elements of P_HT: %d\n", P_HT_coo->N);
	/*for (i = 0; i < P_HT_coo->N; i++) {
		printf("%f ", P_HT_coo->v[i]);
	}
	printf("\n");*/

	/* Export P_HT to .csv file in CSR format */
	printf("Write P_HT to .csv file in CSR format\n");
	sparse_rcs *P_HT;
	P_HT = sparse_rcs_create(P_HT_coo->m, P_HT_coo->n, P_HT_coo->N);
	matrix_trans((const sparse_coo *) P_HT_coo, P_HT);
	FILE *fid;
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
	sparse_coo_free(H);
	sparse_coo_free(C_e_2);
	sparse_coo_free(H_T);
	sparse_coo_free(P_HT_coo);

	sparse_rcs_free(P_HT);

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
	printf("Free all the memory\n");

	return 0;
}