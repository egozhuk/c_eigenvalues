#include "return_codes.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPSILON 1e-7

void free_matrix(double **matrix, int n);
void qr_algorithm(double **A, double **Q, double **R, double **A_next, double *eigenvalues, int n, int max_iterations);
void gram_schmidt(double **A, double **Q, double **R, int n);

int main(int argc, char *argv[])
{
	//	if (argc != 3)
	//	{
	//		printf("WrongArgsError\n");
	//		return ERROR_PARAMETER_INVALID;
	//	}

	FILE *inFile;
	inFile = fopen("/Users/vladislavzukov/Desktop/C:C++/C/lab1/in.txt", "r");
	if (inFile == NULL)
	{
		printf("FileNotFoundError\n");
		fclose(inFile);
		return ERROR_CANNOT_OPEN_FILE;
	}

	int n;
	fscanf(inFile, "%d", &n);

	double **A = (double **)malloc(n * sizeof(double *));

	if (A == NULL)
	{
		printf("MemoryAllocationError\n");
		return ERROR_OUT_OF_MEMORY;
	}

	for (int i = 0; i < n; ++i)
	{
		A[i] = (double *)malloc(n * sizeof(double));

		if (A[i] == NULL)
		{
			printf("MemoryAllocationError\n");
			return ERROR_OUT_OF_MEMORY;
		}

		for (int j = 0; j < n; ++j)
		{
			fscanf(inFile, "%lf ", &A[i][j]);
		}
	}

	if (inFile)
	{
		fclose(inFile);
	}

	//	for (int i = 0; i < n; ++i)
	//	{
	//		for (int j = 0; j < n; ++j)
	//		{
	//			printf("%g ", A[i][j]);
	//		}
	//		printf("\n");
	//	}

	double *eigenvalues = (double *)malloc(n * sizeof(double));
	double **Q = (double **)malloc(n * sizeof(double *));
	double **R = (double **)malloc(n * sizeof(double *));
	double **A_next = (double **)malloc(n * sizeof(double *));

	if (Q == NULL || R == NULL || A_next == NULL || eigenvalues == NULL)
	{
		printf("MemoryAllocationError\n");
		return ERROR_OUT_OF_MEMORY;
	}

	for (int i = 0; i < n; ++i)
	{
		Q[i] = (double *)malloc(n * sizeof(double));
		R[i] = (double *)malloc(n * sizeof(double));
		A_next[i] = (double *)malloc(n * sizeof(double));

		if (A_next[i] == NULL || Q[i] == NULL || R[i] == NULL)
		{
			printf("MemoryAllocationError\n");
			return ERROR_OUT_OF_MEMORY;
		}
	}

	int max_iterations = 10;
	qr_algorithm(A, Q, R, A_next, eigenvalues, n, max_iterations);

	FILE *outFile;
	outFile = fopen("/Users/vladislavzukov/Desktop/C:C++/C/lab1/out.txt", "w");

	if (outFile == NULL)
	{
		printf("FileNotFoundError\n");
		fclose(outFile);
		return ERROR_CANNOT_OPEN_FILE;
	}

	for (int i = 0; i < n; ++i)
	{
		fprintf(outFile, "%f\n", eigenvalues[i]);
	}

	if (outFile)
	{
		fclose(outFile);
	}

	free_matrix(A_next, n);
	free_matrix(Q, n);
	free_matrix(R, n);
	free_matrix(A, n);
	free(eigenvalues);
	return SUCCESS;
}

void free_matrix(double **matrix, int n)
{
	for (int i = 0; i < n; ++i)
	{
		free(matrix[i]);
	}
	free(matrix);
}

//void gram_schmidt(double **A, double **Q, double **R, int n)
//{
//	int i, j, k;
//	for (i = 0; i < n; i++)
//	{
//		for (j = 0; j < n; j++)
//		{
//			R[i][j] = 0;
//		}
//	}
//	for (i = 0; i < n; i++)
//	{
//		double norm = 0;
//		for (j = 0; j < n; j++)
//		{
//			norm += A[j][i] * A[j][i];
//		}
//		norm = sqrt(norm);
//		R[i][i] = norm;
//		for (j = 0; j < n; j++)
//		{
//			Q[j][i] = A[j][i] / norm;
//		}
//		for (j = i + 1; j < n; j++)
//		{
//			double inner_product = 0;
//			for (k = 0; k < n; k++)
//			{
//				inner_product += A[k][j] * Q[k][i];
//			}
//			R[i][j] = inner_product;
//			for (k = 0; k < n; k++)
//			{
//				A[k][j] -= R[i][j] * Q[k][i];
//			}
//		}
//	}
//}

void gram_schmidt(double **A, double **Q, double **R, int n)
{
	int i, j, k;
	double dot_product;
	for (i = 0; i < n; ++i)
	{
		// Set the i-th column of Q to the i-th column of A
		for (j = 0; j < n; ++j)
		{
			Q[j][i] = A[j][i];
		}

		// Subtract the projections of the previous columns from the current column
		for (j = 0; j < i; ++j)
		{
			dot_product = 0;
			for (k = 0; k < n; ++k)
			{
				dot_product += Q[k][j] * A[k][i];
			}
			for (k = 0; k < n; ++k)
			{
				Q[k][i] -= dot_product * Q[k][j];
			}
		}
		// Normalize the current column of Q
		double norm = 0;
		for (j = 0; j < n; ++j)
		{
			norm += Q[j][i] * Q[j][i];
		}
		norm = sqrt(norm);
		for (j = 0; j < n; ++j)
		{
			Q[j][i] /= norm;
		}
		// Calculate the i-th row of R
		for (j = 0; j <= i; ++j)
		{
			dot_product = 0;
			for (k = 0; k < n; ++k)
			{
				dot_product += Q[k][j] * A[k][i];
			}
			R[j][i] = dot_product;
		}
	}
}

void qr_algorithm(double **A, double **Q, double **R, double **A_next, double *eigenvalues, int n, int max_iterations)
{
	for (int iteration = 0; iteration < max_iterations; ++iteration)
	{
		gram_schmidt(A, Q, R, n);

		// A_next = R * Q
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				A_next[i][j] = 0;
				for (int k = 0; k < n; ++k)
				{
					A_next[i][j] += R[i][k] * Q[k][j];
				}
			}
		}

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				A[i][j] = A_next[i][j];
			}
		}

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				printf("%g ", A[i][j]);
			}
			printf("\n");
		}
	}

	for (int i = 0; i < n; ++i)
	{
		eigenvalues[i] = A_next[i][i];
	}
}
