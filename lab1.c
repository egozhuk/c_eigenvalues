#include "return_codes.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void gram_schmidt(int n, double **A, double **Q, double **R)
{
	int i, j, k;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			R[i][j] = 0;
		}
	}

	for (i = 0; i < n; i++)
	{
		double norm = 0;
		for (j = 0; j < n; j++)
		{
			norm += A[j][i] * A[j][i];
		}
		norm = sqrt(norm);
		R[i][i] = norm;

		for (j = 0; j < n; j++)
		{
			Q[j][i] = A[j][i] / norm;
		}

		for (j = i + 1; j < n; j++)
		{
			double inner_product = 0;
			for (k = 0; k < n; k++)
			{
				inner_product += A[k][j] * Q[k][i];
			}
			R[i][j] = inner_product;

			for (k = 0; k < n; k++)
			{
				A[k][j] -= R[i][j] * Q[k][i];
			}
		}
	}
}

void qr_algorithm(double **A, int n, double tol, double **eigenvalues, double **Q, double **R)
{
	const int MAX_ITERATIONS = 8000;
	for (int l = 0; l < MAX_ITERATIONS; l++)
	{
		gram_schmidt(n, A, Q, R);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				double temp = 0;
				for (int k = 0; k < n; k++)
				{
					temp += R[i][k] * Q[k][j];
				}
				A[i][j] = temp;
			}
		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (fabs(A[i][j]) < tol)
			{
				A[i][j] = 0;
			}
		}
	}

	for (int i = 0; i < n; i++)
	{
		if (i != n - 1 && A[i + 1][i] != 0)
		{
			eigenvalues[i][0] = (A[i][i] + A[i + 1][i + 1]) / 2;
			eigenvalues[i + 1][0] = (A[i][i] + A[i + 1][i + 1]) / 2;
			double sum = fabs(A[i + 1][i]) + fabs(A[i][i + 1]);
			eigenvalues[i][1] = (sum) / 2 * -1;
			eigenvalues[i + 1][1] = sum / 2;
			i++;
		}
		else
		{
			eigenvalues[i][0] = A[i][i];
			eigenvalues[i][1] = 0;
		}
	}
}

int main(int argc, char *argv[])
{
	int res = SUCCESS;
	if (argc != 3)
	{
		fprintf(stderr, "Wrong number of program arguments");
		return ERROR_PARAMETER_INVALID;
	}
	int n;
	FILE *in = NULL;
	FILE *out = NULL;
	in = fopen(argv[1], "r");
	if (in == NULL)
	{
		fprintf(stderr, "File %s not found\n", argv[1]);
		return ERROR_CANNOT_OPEN_FILE;
	}
	fscanf(in, "%d", &n);
	const double tol = 1e-3;
	double **A = (double **)malloc(n * sizeof(double *));
	if (A == NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		res = ERROR_OUT_OF_MEMORY;
		goto freeMem;
	}
	for (int i = 0; i < n; i++)
	{
		A[i] = (double *)malloc(n * sizeof(double));
		if (A[i] == NULL)
		{
			fprintf(stderr, "Memory allocation failed\n");
			res = ERROR_OUT_OF_MEMORY;
			goto freeMem;
		}
		for (int j = 0; j < n; j++)
		{
			fscanf(in, "%lf", &A[i][j]);
		}
	}

	double **Q = (double **)malloc(n * sizeof(double *));
	double **R = (double **)malloc(n * sizeof(double *));
	if (Q == NULL || R == NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		res = ERROR_OUT_OF_MEMORY;
		goto freeMem;
	}
	for (int i = 0; i < n; i++)
	{
		Q[i] = (double *)malloc(n * sizeof(double));
		R[i] = (double *)malloc(n * sizeof(double));
		if (Q[i] == NULL || R[i] == NULL)
		{
			fprintf(stderr, "Memory allocation failed\n");
			res = ERROR_OUT_OF_MEMORY;
			goto freeMem;
		}
	}

	double **eigenvalues = (double **)malloc(n * sizeof(double *));
	if (eigenvalues == NULL)
	{
		fprintf(stderr, "Memory allocation failed\n");
		res = ERROR_OUT_OF_MEMORY;
		goto freeMem;
	}
	for (int i = 0; i < n; i++)
	{
		eigenvalues[i] = (double *)malloc(2 * sizeof(double));
		if (eigenvalues[i] == NULL)
		{
			fprintf(stderr, "Memory allocation failed\n");
			res = ERROR_OUT_OF_MEMORY;
			goto freeMem;
		}
	}
	qr_algorithm(A, n, tol, eigenvalues, Q, R);

	out = fopen(argv[2], "w");
	if (out == NULL)
	{
		fprintf(stderr, "File %s not found\n", argv[2]);
		res = ERROR_CANNOT_OPEN_FILE;
		goto freeMem;
	}
	for (int i = 0; i < n; i++)
	{
		if (eigenvalues[i][1] == 0)
		{
			fprintf(out, "%g\n", eigenvalues[i][0]);
		}
		else
		{
			if (eigenvalues[i][1] >= 0)
			{
				fprintf(out, "%g +%gi\n", eigenvalues[i][0], eigenvalues[i][1]);
			}
			else
			{
				fprintf(out, "%g %gi\n", eigenvalues[i][0], eigenvalues[i][1]);
			}
		}
	}

freeMem:
	if (in)
	{
		fclose(in);
	}
	if (out)
	{
		fclose(out);
	}

	if (A != NULL)
	{
		for (int i = 0; i < n; i++)
		{
			if (A[i] != NULL)
				free(A[i]);
		}
		free(A);
	}

	if (Q != NULL)
	{
		for (int i = 0; i < n; i++)
		{
			if (Q[i] != NULL)
				free(Q[i]);
		}
		free(Q);
	}

	if (R != NULL)
	{
		for (int i = 0; i < n; i++)
		{
			if (R[i] != NULL)
				free(R[i]);
		}
		free(R);
	}

	if (eigenvalues != NULL)
	{
		for (int i = 0; i < n; i++)
		{
			if (eigenvalues[i] != NULL)
				free(eigenvalues[i]);
		}
		free(eigenvalues);
	}

	return res;
}


//		printf("[");
//		for (int i = 0; i < n; ++i)
//		{
//			if (i + 1 == n)
//			{
//				printf("%f, %f, %f, %f]", Q[i][0], Q[i][1], Q[i][2], Q[i][3]);
//				break;
//			}
//			printf("%f, %f, %f, %f; ", Q[i][0], Q[i][1], Q[i][2], Q[i][3]);
//		}
//		printf("\n");
//		printf("[");
//		for (int i = 0; i < n; ++i)
//		{
//			if (i + 1 == n)
//			{
//				printf("%f, %f, %f, %f]", R[i][0], R[i][1], R[i][2], R[i][3]);
//				break;
//			}
//			printf("%f, %f, %f, %f; ", R[i][0], R[i][1], R[i][2], R[i][3]);
//		}
//		printf("\n");
//		return;