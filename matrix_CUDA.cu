
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <device_functions.h>
#include <stdio.h>
#include <stdlib.h>

#include "type.h"

#define TILE_WIDTH 32
__global__ void denseMatrixTranspose(float *A, float *A_T, int numARows, int numAColumns) {
	// The block and thread indices
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Identify the row and column of the P element
	int Row = by * TILE_WIDTH + ty;
	int Col = bx * TILE_WIDTH + tx;

	if (Row < numARows && Col < numAColumns) {
		A_T[Col*numARows + Row] = A[Row*numAColumns + Col];
	}
}

extern "C" void denseMatrixTranspose_wrapper(float **A, float **A_T, int numARows, int numAColumns) {
	float *d_A, *d_A_T;
	int i;

	/* Allocate device memory and copy host memory to device memory */
	cudaMalloc((void **)&d_A, numARows*numAColumns * sizeof(float));
	cudaMalloc((void **)&d_A_T, numARows*numAColumns * sizeof(float));

	/* Copy from host to device */
	for (i = 0; i < numARows; i++) {
		cudaMemcpy(d_A + (i*numAColumns), A[i], numAColumns * sizeof(float), cudaMemcpyHostToDevice);
	}

	/* Initialize the grid and block dimensions */
	dim3 DimGrid(ceil(1.0*numAColumns / TILE_WIDTH), ceil(1.0*numARows / TILE_WIDTH), 1);
	dim3 DimBlock(TILE_WIDTH, TILE_WIDTH, 1);

	/* Launch the GPU kernel */
	denseMatrixTranspose<<< DimGrid, DimBlock >>>(d_A, d_A_T, numARows, numAColumns);
	cudaDeviceSynchronize();
	/* Copyt the GPU memory back to the CPU */
	for (i = 0; i < numAColumns; i++) {
		cudaMemcpy(A_T[i], d_A_T + (i*numARows), numARows * sizeof(float), cudaMemcpyDeviceToHost);
	}

	/* Free the GPU memory */
	cudaFree(d_A);
	cudaFree(d_A_T);
}

__global__ void denseMatrixMultiply(float *A, float *B, float *C, int numARows, int numAColumns, int numBRows, int numBColumns) {
	__shared__ float subTileA[TILE_WIDTH][TILE_WIDTH];
	__shared__ float subTileB[TILE_WIDTH][TILE_WIDTH];

	// The block and thread indices
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Identify the row and column of the P element
	int Row = by * TILE_WIDTH + ty;
	int Col = bx * TILE_WIDTH + tx;
	float Pvalue = 0;

	// Loop over the M and N tiles required to compute the P element
	for (int m = 0; m < (numAColumns - 1) / TILE_WIDTH + 1; m++) {
		// Collaborative loading of M and N tiles into shared memory
		if (Row < numARows && m*TILE_WIDTH + tx < numAColumns) {
			subTileA[ty][tx] = A[Row*numAColumns + m*TILE_WIDTH + tx];
		}
		else {
			subTileA[ty][tx] = 0;
		}
		if (m*TILE_WIDTH + ty < numBRows && Col < numBColumns) {
			subTileB[ty][tx] = B[(m*TILE_WIDTH + ty)*numBColumns + Col];
		}
		else {
			subTileB[ty][tx] = 0;
		}

		__syncthreads();
		if (Row < numARows && Col < numBColumns) {
			for (int k = 0; k < TILE_WIDTH; k++) {
				Pvalue += subTileA[ty][k] * subTileB[k][tx];
			}
		}

		__syncthreads();
	}

	if (Row < numARows && Col < numBColumns)
		C[Row*numBColumns + Col] = Pvalue;
}

extern "C" void denseMatrixMultiply_wrapper(float **A, float **B, float **C, int numARows, int numAColumns, int numBRows, int numBColumns) {
	float *d_A, *d_B, *d_C;
	int i;

	/* Allocate device memory and copy host memory to device memory */
	cudaMalloc((void **)&d_A, numARows*numAColumns * sizeof(float));
	cudaMalloc((void **)&d_B, numBRows*numBColumns * sizeof(float));
	cudaMalloc((void **)&d_C, numARows*numBColumns * sizeof(float));

	/* Copy from host to device */
	for (i = 0; i < numARows; i++) {
		cudaMemcpy(d_A + (i*numAColumns), A[i], numAColumns * sizeof(float), cudaMemcpyHostToDevice);
	}
	for (i = 0; i < numBRows; i++) {
		cudaMemcpy(d_B + (i*numBColumns), B[i], numBColumns * sizeof(float), cudaMemcpyHostToDevice);
	}

	/* Initialize the grid and block dimensions */
	dim3 DimGrid(ceil(1.0*numBColumns / TILE_WIDTH), ceil(1.0*numARows / TILE_WIDTH), 1);
	dim3 DimBlock(TILE_WIDTH, TILE_WIDTH, 1);

	/* Launch the GPU kernel */
	denseMatrixMultiply<<< DimGrid, DimBlock >>>(d_A, d_B, d_C, numARows, numAColumns, numBRows, numBColumns);
	cudaDeviceSynchronize();
	/* Copyt the GPU memory back to the CPU */
	for (i = 0; i < numARows; i++) {
		cudaMemcpy(C[i], d_C + (i*numBColumns), numBColumns * sizeof(float), cudaMemcpyDeviceToHost);
	}

	/* Free the GPU memory */
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
}

__global__ void sparseVector(const int num_rows,
	const int *ptr,
	const int * indices,
	const float * data,
	const float * x,
	float * y,
	int Div)
{

	//__shared__ float vals[num_rows * 32];

	//int tid = blockDim.x*blockIdx.x + threadIdx.x;
	//int warp_id = tid / 32;
	//int lane = tid & (32 - 1);

	//int row = warp_id;

	//if (row < num_rows) {
	//	int row_start = ptr[row];
	//	int row_end = ptr[row + 1];

	//	vals[threadIdx.x] = 0;
	//	for (int jj = row_start + lane; jj < row_end; jj += 32)
	//		vals[threadIdx.x] += data[jj] * x[indices[jj]];

	//	if (lane < 16) vals[threadIdx.x] += vals[threadIdx.x + 16];
	//	if (lane < 8) vals[threadIdx.x] += vals[threadIdx.x + 8];
	//	if (lane < 4) vals[threadIdx.x] += vals[threadIdx.x + 4];
	//	if (lane < 2) vals[threadIdx.x] += vals[threadIdx.x + 2];
	//	if (lane < 1) vals[threadIdx.x] += vals[threadIdx.x + 1];

	//	if (lane == 0) {
	//		y[row] += vals[threadIdx.x]/Div;
	//	}
	//}


	   //printf("%f\n", data[0]);
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	
	if(row < num_rows){
	    float dot = 0;
	    int row_start = ptr[row];
	    int row_end = ptr[row+1];
	
	    for (int jj = row_start; jj < row_end; jj++)
	    dot+= data[jj] * x[indices[jj]];
	
	    y[row] += dot / Div;
	   // printf("%f\n", y[row]);
	}
}

extern "C" void sparseVector_wrapper(sparse_rcs *A, float **B, float **C, int FR_RC, int Div) {
	float *d_v, *d_B, *d_C;
	int *d_j, *d_r;
	int i;
	/* Allocate device memory and copy host memory to device memory */
	cudaMalloc((void **)&d_v, A->N * sizeof(float));
	cudaMalloc((void **)&d_j, A->N * sizeof(int));
	cudaMalloc((void **)&d_r, A->N * sizeof(int));
	cudaMalloc((void **)&d_B, (FR_RC*FR_RC) * sizeof(float));
	cudaMalloc((void **)&d_C, (FR_RC*A->m) * sizeof(float));

	/* Copy from host to device */
	cudaMemcpy(d_v, A->v, A->N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, A->r, A->N * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_j, A->j, A->N * sizeof(int), cudaMemcpyHostToDevice);
	for (i = 0; i < FR_RC; i++) {
		cudaMemcpy(d_B + (i*FR_RC), B[i], FR_RC * sizeof(float), cudaMemcpyHostToDevice);
	}

	/* Initialize the grid and block dimensions */
	dim3 DimGrid(ceil(1.0*FR_RC / TILE_WIDTH), 1, 1);
	dim3 DimBlock(TILE_WIDTH, 1, 1);

	/* Launch the GPU kernel */
	for (i = 0; i < FR_RC; i++) {
		sparseVector <<< DimGrid, DimBlock >>> (A->m, d_r, d_j, d_v, d_B + (i*FR_RC), d_C + (i*FR_RC), Div);
		cudaDeviceSynchronize();
	}
	/* Copyt the GPU memory back to the CPU */
	for (i = 0; i < FR_RC; i++) {
		cudaMemcpy(C[i], d_C + (i*FR_RC), FR_RC * sizeof(float), cudaMemcpyDeviceToHost);
	}

	/* Free the GPU memory */
	cudaFree(d_v);
	cudaFree(d_r);
	cudaFree(d_j);
	cudaFree(d_B);
	cudaFree(d_C);
}

__global__ void cooToFullRank(float *v, int *i, int *j, float *Out, int len, int numCol) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < len) {
		Out[i[idx] * numCol + j[idx]] = v[idx];
	}
}

extern "C" void cooToFullRank_wrapper(sparse_coo *A_coo, float **A_FR) {
	float *d_v, *d_out;
	int *d_i, *d_j;
	int i;
	/* Allocate device memory and copy host memory to device memory */
	cudaMalloc((void **)&d_v, A_coo->N * sizeof(float));
	cudaMalloc((void **)&d_i, A_coo->N * sizeof(int));
	cudaMalloc((void **)&d_j, A_coo->N * sizeof(int));
	cudaMalloc((void **)&d_out, (A_coo->m*A_coo->n) * sizeof(float));

	/* Copy from host to device */
	cudaMemcpy(d_v, A_coo->v, A_coo->N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_i, A_coo->i, A_coo->N * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_j, A_coo->j, A_coo->N * sizeof(int), cudaMemcpyHostToDevice);

	/* Initialize the grid and block dimensions */
	dim3 DimGrid(ceil(1.0*A_coo->N / TILE_WIDTH), 1, 1);
	dim3 DimBlock(TILE_WIDTH, 1, 1);

	/* Launch the GPU kernel */
	cooToFullRank <<< DimGrid, DimBlock >>> (d_v, d_i, d_j, d_out, A_coo->N, A_coo->n);
	cudaDeviceSynchronize();

	/* Copyt the GPU memory back to the CPU */
	for (i = 0; i < A_coo->m; i++) {
		cudaMemcpy(A_FR[i], d_out + (i*A_coo->n), A_coo->n * sizeof(float), cudaMemcpyDeviceToHost);
	}

	/* Free the GPU memory */
	cudaFree(d_v);
	cudaFree(d_i);
	cudaFree(d_j);
	cudaFree(d_out);
}