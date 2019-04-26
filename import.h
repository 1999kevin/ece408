#ifndef IMPORT_H
#define IMPORT_H

#include <stdio.h>
#include "type.h"

sparse_rcs *sparse_rcs_create(int m, int n, int N);

void sparse_rcs_free(sparse_rcs *A);

sparse_coo *sparse_coo_create(int m, int n, int N);

void sparse_coo_free(sparse_coo *A);

sparse_coo *sparse_coo_import(const char *filename);

full_r *full_r_create(int m, int n);

void full_r_free(full_r *A);

full_r *full_r_import(const char *filename);

toeplitz *toeplitz_create(int N, int n);

void toeplitz_free(toeplitz *A);

toeplitz *toeplitz_import(const char *filename);

#endif
