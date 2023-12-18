// #include <complex.h>
// #include <math.h>
// #include <stdio.h>
//
// #define N 2
//
// typedef complex double Matrix[N][N];
//
// void householder(Matrix A, Matrix Q, Matrix R) {
//	for (int k = 0; k < N; k++) {
//		complex double x[N - k];
//		complex double alpha = 0;
//		for (int i = 0; i < N - k; i++) {
//			x[i] = A[k + i][k];
//			alpha += x[i] * conj(x[i]);
//		}
//		alpha = sqrt(alpha);
//		if (cabs(x[0]) > 1e-9) {
//			x[0] -= copysign(alpha, creal(x[0])) + I * cimag(x[0]);
//		} else {
//			x[0] -= alpha;
//		}
//		complex double beta = 0;
//		for (int i = 0; i < N - k; i++) {
//			beta += x[i] * conj(x[i]);
//		}
//		beta = sqrt(beta);
//		for (int i = 0; i < N - k; i++) {
//			x[i] /= beta;
//		}
//		for (int i = 0; i < N - k; i++) {
//			for (int j = 0; j < N - k; j++) {
//				complex double v = x[i] * conj(x[j]);
//				A[k + i][k + j] -= 2 * v;
//			}
//		}
//	}
//
//	for (int i = 0; i < N; i++) {
//		for (int j = 0; j < N; j++) {
//			if (i > j) {
//				Q[i][j] = A[i][j];
//				R[i][j] = 0;
//			} else {
//				Q[i][j] = 0;
//				R[i][j] = A[i][j];
//			}
//		}
//		Q[i][i] += 1;
//	}
// }
//
// void schur(Matrix A, Matrix Q, Matrix R) {
//	Matrix temp;
//	while (true) {
//		bool converged = true;
//		for (int i = 0; i < N - 1; i++) {
//			if (cabs(A[i + 1][i]) > 1e-9) {
//				converged = false;
//				break;
//			}
//		}
//		if (converged) {
//			break;
//		}
//		householder(A, Q, R);
//		// Multiply A = Q^H * A
//		for (int i = 0; i < N; i++) {
//			for (int j = 0; j < N; j++) {
//				temp[i][j] = 0;
//				for (int k = 0; k < N; k++) {
//					temp[i][j] += conj(Q[k][i]) * A[k][j];
//				}
//			}
//		}
//		// Multiply A = A * R
//		for (int i = 0; i < N; i++) {
//			for (int j = 0; j < N; j++) {
//				A[i][j] = 0;
//				for (int k = 0; k < N; k++) {
//					A[i][j] += temp[i][k] * R[k][j];
//				}
//			}
//		}
//	}
// }
//
// int main() {
//	Matrix A = {
//		{ 0 + 0*I, 1 + 0*I },
//		{ -1 - 2*I, 0 + 0*I },
//	};
//
//	Matrix Q, R;
//	schur(A, Q, R);
//
//	printf("Eigenvalues:\n");
//	for (int i = 0; i < N; i++) {
//		printf("%g + %gi\n", creal(A[i][i]), cimag(A[i][i]));
//	}
//
//	return 0;
// }

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
//
// void hessenberg_reduction(double **A, int n);
// void qr_decomposition(double **A, double **Q, double **R, int n);
// void qr_algorithm(double **A, double *eigenvalues, int n, int max_iterations, double epsilon);
//
// int main() {
//	int n = 2;
//	double **A = (double **)malloc(n * sizeof(double *));
//	for (int i = 0; i < n; i++) {
//		A[i] = (double *)malloc(n * sizeof(double));
//	}
//
//	A[0][0] = 0; A[0][1] = 1.0;
//	A[1][0] = -1.0; A[1][1] = 0;
//
//	double *eigenvalues = (double *)malloc(n * sizeof(double));
//	int max_iterations = 1000;
//	double epsilon = 1e-8;
//	qr_algorithm(A, eigenvalues, n, max_iterations, epsilon);
//
//	printf("Eigenvalues:\n");
//	for (int i = 0; i < n; i++) {
//		printf("%g\n", eigenvalues[i]);
//	}
//
//	for (int i = 0; i < n; i++) {
//		free(A[i]);
//	}
//	free(A);
//	free(eigenvalues);
//
//	return 0;
// }
//
// void hessenberg_reduction(double **A, int n) {
//	for (int k = 0; k < n - 2; k++) {
//		double x = 0;
//		int m = k + 1;
//		for (int i = k + 1; i < n; i++) {
//			if (fabs(A[i][k]) > x) {
//				x = fabs(A[i][k]);
//				m = i;
//			}
//		}
//
//		if (m != k + 1) {
//			for (int j = k; j < n; j++) {
//				double temp = A[m][j];
//				A[m][j] = A[k+1][j];
//				A[k+1][j] = temp;
//			}
//
//			for (int i = 0; i < n; i++) {
//				double temp = A[i][m];
//				A[i][m] = A[i][k+1];
//				A[i][k+1] = temp;
//			}
//		}
//
//		if (fabs(A[k+1][k]) > 1e-9) {
//			double t = sqrt(2 * (A[k+1][k] * A[k+1][k]));
//			double r = A[k+1][k] > 0 ? -A[k+1][k] - t : -A[k+1][k] + t;
//			double c = sqrt(0.5 * (1 - A[k+1][k] / r));
//			double s = -A[k+2][k] / (2 * c * r);
//
//			for (int j = k + 1; j < n; j++) {
//				double temp = c * A[k+1][j] - s * A[k+2][j];
//				A[k+2][j] = s * A[k+1][j] + c * A[k+2][j];
//				A[k+1][j] = temp;
//			}
//
//			for (int i = 0; i < n; i++) {
//				double temp = c * A[i][k+1] - s * A[i][k+2];
//				A[i][k+2] = s * A[i][k+1] + c * A[i][k+2];
//				A[i][k+1] = temp;
//			}
//		}
//	}
// }
//
// void qr_decomposition(double **A, double **Q, double **R, int n) {
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			R[i][j] = A[i][j];
//		}
//	}
//
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			Q[i][j] = (i == j) ? 1.0 : 0.0;
//		}
//	}
//
//	for (int k = 0; k < n - 1; k++) {
//		double x = 0;
//		for (int i = k + 1; i < n; i++) {
//			x += R[i][k] * R[i][k];
//		}
//		double r = sqrt(R[k][k] * R[k][k] + x);
//		double c = R[k][k] / r;
//		double s = -sqrt(x) / r;
//
//		for (int i = k; i < n; i++) {
//			double temp = c * R[k][i] - s * R[k+1][i];
//			R[k+1][i] = s * R[k][i] + c * R[k+1][i];
//			R[k][i] = temp;
//		}
//
//		for (int i = 0; i < n; i++) {
//			double temp = c * Q[i][k] - s * Q[i][k+1];
//			Q[i][k+1] = s * Q[i][k] + c * Q[i][k+1];
//			Q[i][k] = temp;
//		}
//	}
// }
//
// void qr_algorithm(double **A, double *eigenvalues, int n, int max_iterations, double epsilon) {
//	hessenberg_reduction(A, n);
//	double **Q = (double **)malloc(n * sizeof(double *));
//	double **R = (double **)malloc(n * sizeof(double *));
//	for (int i = 0; i < n; i++) {
//		Q[i] = (double *)malloc(n * sizeof(double));
//		R[i] = (double *)malloc(n * sizeof(double));
//	}
//
//	for (int iter = 0; iter < max_iterations; iter++) {
//		qr_decomposition(A, Q, R, n);
//
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < n; j++) {
//				A[i][j] = 0;
//				for (int k = 0; k < n; k++) {
//					A[i][j] += R[i][k] * Q[k][j];
//				}
//			}
//		}
//
//		double max_off_diagonal = 0;
//		for (int i = 0; i < n - 1; i++) {
//			max_off_diagonal = fmax(max_off_diagonal, fabs(A[i][i + 1]));
//		}
//
//		if (max_off_diagonal < epsilon) {
//			break;
//		}
//	}
//
//	for (int i = 0; i < n; i++) {
//		eigenvalues[i] = A[i][i];
//	}
//
//	for (int i = 0; i < n; i++) {
//		free(Q[i]);
//		free(R[i]);
//	}
//	free(Q);
//	free(R);
// }

//#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
//
//typedef struct
//{
//	double real;
//	double imag;
//} Complex;
//
//Complex complex_add(Complex a, Complex b);
//Complex complex_sub(Complex a, Complex b);
//Complex complex_mul(Complex a, Complex b);
//Complex complex_div(Complex a, Complex b);
//double complex_abs(Complex a);
//
//void householder(Complex **A, Complex **Q, Complex **R, int n);
//void qr_algorithm(Complex **A, int n, double tol, Complex *eigenvalues);
//
//int main(int argc, char **argv)
//{
//	int n;
//	FILE *f;
//	f = fopen(argv[1], "r");
//	fscanf(f, "%d", &n);
//	double tol = 1e-12;
//	Complex **A = malloc(sizeof(Complex *) * n);
//	for (int i = 0; i < n; i++)
//	{
//		A[i] = malloc(sizeof(Complex) * n);
//		for (int j = 0; j < n; j++)
//		{
//			fscanf(f, "%f", &A[i][j].real);
//			A[i][j].imag = 0;
//		}
//	}
//	fclose(f);
//	Complex *eigenvalues = malloc(n * sizeof(Complex));
//
//	qr_algorithm(A, n, tol, eigenvalues);
//
//	f = fopen(argv[2], "w");
//	for (int i = 0; i < n; i++)
//	{
//		fprintf(f, "%g +%gi\n", eigenvalues[i].real, eigenvalues[i].imag);
//	}
//	fclose(f);
//	for (int i = 0; i < n; ++i)
//	{
//		free(A[i]);
//	}
//	free(A);
//	free(eigenvalues);
//	return 0;
//}
//
//Complex complex_add(Complex a, Complex b)
//{
//	Complex result;
//	result.real = a.real + b.real;
//	result.imag = a.imag + b.imag;
//	return result;
//}
//
//Complex complex_sub(Complex a, Complex b)
//{
//	Complex result;
//	result.real = a.real - b.real;
//	result.imag = a.imag - b.imag;
//	return result;
//}
//
//Complex complex_mul(Complex a, Complex b)
//{
//	Complex result;
//	result.real = a.real * b.real - a.imag * b.imag;
//	result.imag = a.real * b.imag + a.imag * b.real;
//	return result;
//}
//
//Complex complex_div(Complex a, Complex b)
//{
//	Complex result;
//	double denominator = b.real * b.real + b.imag * b.imag;
//	result.real = (a.real * b.real + a.imag * b.imag) / denominator;
//	result.imag = (a.imag * b.real - a.real * b.imag) / denominator;
//	return result;
//}
//
//Complex complex_sqrt(Complex a) {
//	Complex b = {0, 0};
//	if (a.imag != 0) {
//		double mag = sqrt(a.real * a.real + a.imag * a.imag);
//		double angle = atan2(a.imag, a.real);
//
//		b.real = sqrt(mag) * cos(angle / 2);
//		b.imag = sqrt(mag) * sin(angle / 2);
//	} else if (a.real >= 0) {
//		b.real = sqrt(a.real);
//	} else {
//		b.imag = sqrt(-1 * a.real);
//	}
//	return b;
//}
//
//Complex complex_copysign(Complex x, Complex y)
//{
//	double sign_re = y.real < 0 ? -1.0 : 1.0;
//	double sign_im = y.imag < 0 ? -1.0 : 1.0;
//	return (Complex){ sign_re * x.real, sign_im * x.imag };
//}
//
//double complex_abs(Complex a)
//{
//	return sqrt(a.real * a.real + a.imag * a.imag);
//}
//
//void householder(Complex **A, Complex **Q, Complex **R, int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			R[i][j] = A[i][j];
//			if (i == j)
//			{
//				Q[i][j].real = 1;
//				Q[i][j].imag = 0;
//			}
//			else
//			{
//				Q[i][j].real = 0;
//				Q[i][j].imag = 0;
//			}
//		}
//	}
//
//	for (int k = 0; k < n - 1; k++)
//	{
//		Complex x_norm = { 0.0, 0.0 };
//		for (int i = k; i < n; i++)
//		{
//			x_norm = complex_add(x_norm, complex_mul(R[i][k], R[i][k]));
//		}
//		x_norm = complex_sqrt(x_norm);
//
//		Complex v[n];
//		for (int i = 0; i < n; i++)
//		{
//			if (i < k)
//			{
//				v[i].real = 0;
//				v[i].imag = 0;
//			}
//			else if (i == k)
//			{
//				v[i] = complex_add(R[i][k], complex_copysign(x_norm, R[i][k]));
//			}
//			else
//			{
//				v[i] = R[i][k];
//			}
//		}
//
//		Complex v_norm = { 0.0, 0.0 };
//		for (int i = k; i < n; i++)
//		{
//			v_norm = complex_add(v_norm, complex_mul(v[i], v[i]));
//		}
//		v_norm = complex_sqrt(v_norm);
//
//		for (int i = k; i < n; i++)
//		{
//			v[i] = complex_div(v[i], v_norm);
//		}
//
//		for (int i = k; i < n; i++)
//		{
//			Complex temp = { 0.0, 0.0 };
//			for (int j = k; j < n; j++)
//			{
//				temp = complex_add(temp, complex_mul(R[j][i], v[j]));
//			}
//			for (int j = k; j < n; j++)
//			{
//				Complex mul = { 2.0, 0.0 };
//				R[j][i] = complex_sub(R[j][i], complex_mul(mul, complex_mul(temp, v[j])));
//			}
//		}
//
//		for (int i = 0; i < n; i++)
//		{
//			Complex temp = { 0.0, 0.0 };
//			for (int j = k; j < n; j++)
//			{
//				temp = complex_add(temp, complex_mul(Q[i][j], v[j]));
//			}
//			for (int j = k; j < n; j++)
//			{
//				Complex mul = { 2.0, 0.0 };
//				Q[i][j] = complex_sub(Q[i][j], complex_mul(mul, complex_mul(temp, v[j])));
//			}
//		}
//	}
//}
//
//void qr_algorithm(Complex **A, int n, double tol, Complex *eigenvalues)
//{
//	Complex **Q = malloc(n * sizeof(Complex *));
//	Complex **R = malloc(n * sizeof(Complex *));
//	for (int i = 0; i < n; i++)
//	{
//		Q[i] = malloc(n * sizeof(Complex));
//		R[i] = malloc(n * sizeof(Complex));
//	}
//
//	double err = tol + 1;
//	int iter = 0;
//	int max_iter = 1000;
//
//	while (err > tol && iter < max_iter)
//	{
//		err = 0;
//		// Compute the Rayleigh quotient shift
//		Complex shift = A[n - 1][n - 1];
//
//		// Shift the matrix
//		for (int i = 0; i < n; i++)
//		{
//			A[i][i] = complex_sub(A[i][i], shift);
//		}
//
//		householder(A, Q, R, n);
//
//		// Unshift the matrix and update A
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				Complex temp = { 0.0, 0.0 };
//				for (int k = 0; k < n; k++)
//				{
//					temp = complex_add(temp, complex_mul(R[i][k], Q[k][j]));
//				}
//				A[i][j] = temp;
//				if (i != j && complex_abs(A[i][j]) > err)
//				{
//					err = complex_abs(A[i][j]);
//				}
//			}
//			A[i][i] = complex_add(A[i][i], shift);
//		}
//
//		iter++;
//	}
//
//	if (iter == max_iter)
//	{
//		fprintf(stderr, "QR algorithm did not converge within the specified number of iterations.\n");
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		eigenvalues[i] = A[i][i];
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		free(Q[i]);
//		free(R[i]);
//	}
//	free(Q);
//	free(R);
//}

// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
//
// #define MAX_ITER 1000
// #define EPSILON 1e-12
//
// void qr_decomposition(double **A, double **Q, double **R, int n)
//{
//	for (int k = 0; k < n; k++)
//	{
//		double norm = 0;
//		for (int i = k; i < n; i++)
//		{
//			norm += A[i][k] * A[i][k];
//		}
//		norm = sqrt(norm);
//
//		if (A[k][k] < 0)
//		{
//			norm = -norm;
//		}
//
//		double u = A[k][k] + norm;
//		double gamma = 1 / (u * norm);
//
//		double *w = (double *)calloc(n, sizeof(double));
//		w[k] = u;
//
//		for (int i = k + 1; i < n; i++)
//		{
//			w[i] = A[i][k];
//		}
//
//		for (int j = k; j < n; j++)
//		{
//			double t = 0;
//			for (int i = k; i < n; i++)
//			{
//				t += w[i] * A[i][j];
//			}
//			for (int i = k; i < n; i++)
//			{
//				A[i][j] -= gamma * w[i] * t;
//			}
//		}
//
//		for (int j = 0; j < n; j++)
//		{
//			double t = 0;
//			for (int i = k; i < n; i++)
//			{
//				t += w[i] * Q[i][j];
//			}
//			for (int i = k; i < n; i++)
//			{
//				Q[i][j] -= gamma * w[i] * t;
//			}
//		}
//
//		free(w);
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			R[i][j] = A[i][j];
//		}
//	}
// }
//
// int shifted_qr(double **A, double *eigenvalues, int n)
//{
//	int iter = 0;
//	while (iter < MAX_ITER)
//	{
//		double mu = A[n - 1][n - 1];
//
//		for (int i = 0; i < n; i++)
//		{
//			A[i][i] -= mu;
//		}
//
//		double **Q = (double **)malloc(n * sizeof(double *));
//		double **R = (double **)malloc(n * sizeof(double *));
//		for (int i = 0; i < n; i++)
//		{
//			Q[i] = (double *)malloc(n * sizeof(double));
//			R[i] = (double *)malloc(n * sizeof(double));
//			Q[i][i] = 1.0;
//		}
//
//		qr_decomposition(A, Q, R, n);
//
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				A[i][j] = 0;
//				for (int k = 0; k < n; k++)
//				{
//					A[i][j] += R[i][k] * Q[k][j];
//				}
//			}
//		}
//
//		for (int i = 0; i < n; i++)
//		{
//			A[i][i] += mu;
//		}
//
//		int converged = 1;
//		for (int i = 1; i < n; i++)
//		{
//			if (fabs(A[i][i - 1]) > EPSILON)
//			{
//				converged = 0;
//				break;
//			}
//		}
//
//		if (converged)
//		{
//			break;
//		}
//
//		iter++;
//		for (int i = 0; i < n; i++)
//		{
//			free(Q[i]);
//			free(R[i]);
//		}
//		free(Q);
//		free(R);
//	}
//
//	if (iter == MAX_ITER)
//	{
//		return -1;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		eigenvalues[i] = A[i][i];
//	}
//
//	return 0;
// }
//
// int main(int argc, char **argv)
//{
//	int n;
//	FILE *f;
//	f = fopen(argv[1], "r");
//	fscanf(f, "%d", &n);
//	double **A = (double **)malloc(n * sizeof(double *));
//	for (int i = 0; i < n; i++)
//	{
//		A[i] = malloc(n * sizeof(double));
//		for (int j = 0; j < n; j++)
//		{
//			fscanf(f, "%lf", &A[i][j]);
//		}
//	}
//	fclose(f);
//	double *eigenvalues = malloc(n * sizeof(double));
//
//	shifted_qr(A, eigenvalues, n);
//
//	f = fopen(argv[2], "w");
//	for (int i = 0; i < n; i++)
//	{
//		fprintf(f, "%g\n", eigenvalues[i]);
//	}
//	fclose(f);
//	for (int i = 0; i < n; i++)
//	{
//		free(A[i]);
//	}
//	free(A);
//	free(eigenvalues);
//	return 0;
// }

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void householder(double **A, double **Q, double **R, int n);
void qr_algorithm(double **A, int n, double tol, double *eigenvalues);

int main(int argc, char **argv)
{
	int n;
	FILE *f;
	f = fopen(argv[1], "r");
	fscanf(f, "%d", &n);
	double tol = 1e-12;
	double A[n][n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			fscanf(f, "%lf", &A[i][j]);
		}
	}
	fclose(f);
	double *eigenvalues = malloc(n * sizeof(double));
	double **A_ptr = malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		A_ptr[i] = A[i];
	}

	qr_algorithm(A_ptr, n, tol, eigenvalues);

	f = fopen(argv[2], "w");
	for (int i = 0; i < n; i++)
	{
		fprintf(f, "%g\n", eigenvalues[i]);
	}
	fclose(f);
	free(A_ptr);
	free(eigenvalues);
	return 0;
}

void householder(double **A, double **Q, double **R, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			R[i][j] = A[i][j];
			if (i == j)
			{
				Q[i][j] = 1;
			}
			else
			{
				Q[i][j] = 0;
			}
		}
	}

	for (int k = 0; k < n - 1; k++)
	{
		double x_norm = 0;
		for (int i = k; i < n; i++)
		{
			x_norm += R[i][k] * R[i][k];
		}
		x_norm = sqrt(x_norm);

		double v[n];
		for (int i = 0; i < n; i++)
		{
			if (i < k)
			{
				v[i] = 0;
			}
			else if (i == k)
			{
				v[i] = R[i][k] + copysign(x_norm, R[i][k]);
			}
			else
			{
				v[i] = R[i][k];
			}
		}

		double v_norm = 0;
		for (int i = k; i < n; i++)
		{
			v_norm += v[i] * v[i];
		}
		v_norm = sqrt(v_norm);

		for (int i = k; i < n; i++)
		{
			v[i] /= v_norm;
		}

		for (int i = k; i < n; i++)
		{
			double temp = 0;
			for (int j = k; j < n; j++)
			{
				temp += R[j][i] * v[j];
			}
			for (int j = k; j < n; j++)
			{
				R[j][i] -= 2 * temp * v[j];
			}
		}

		for (int i = 0; i < n; i++)
		{
			double temp = 0;
			for (int j = k; j < n; j++)
			{
				temp += Q[i][j] * v[j];
			}
			for (int j = k; j < n; j++)
			{
				Q[i][j] -= 2 * temp * v[j];
			}
		}
	}
}

void qr_algorithm(double **A, int n, double tol, double *eigenvalues)
{
	double **Q = malloc(n * sizeof(double *));
	double **R = malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		Q[i] = malloc(n * sizeof(double));
		R[i] = malloc(n * sizeof(double));
	}

	double err = tol + 1;

	while (err > tol)
	{
		err = 0;
		householder(A, Q, R, n);
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
				if (i != j && fabs(A[i][j]) > err)
				{
					err = fabs(A[i][j]);
				}
			}
		}
	}

	for (int i = 0; i < n; i++)
	{
		eigenvalues[i] = A[i][i];
	}

	// Free allocated memory
	for (int i = 0; i < n; i++)
	{
		free(Q[i]);
		free(R[i]);
	}
	free(Q);
	free(R);
}











//void hessenberg(double **A, int n)
//{
//	for (int i = 0; i < n - 2; i++)
//	{
//		double *x = (double *)malloc((n - i - 1) * sizeof(double));
//		double w = 0;
//		for (int j = 0; j < (n - i - 1); j++)
//		{
//			x[j] = A[j + 1 + i][i];
//			w += x[j] * x[j];
//		}
//		w = sqrt(w);
//		if (x[0] >= 0)
//		{
//			w *= -1;
//		}
//		double *v = (double *)malloc((n - i - 1) * sizeof(double));
//		v[0] = w - x[0];
//		for (int j = 1; j < (n - i - 1); j++)
//		{
//			v[j] = x[j] * -1;
//		}
//		double **p = (double **)malloc((n - i - 1) * sizeof(double *));
//		double d = 0;
//		for (int j = 0; j < (n - i - 1); j++)
//		{
//			d += v[j] * v[j];
//		}
//		for (int j = 0; j < (n - i - 1); j++)
//		{
//			p[j] = (double *)malloc((n - i - 1) * sizeof(double));
//			for (int k = 0; k < (n - i - 1); k++)
//			{
//				if (d != 0)
//				{
//					p[j][k] = v[k] * v[j] / d;
//				}
//			}
//		}
//		double **h = (double **)malloc((n) * sizeof(double *));
//		for (int j = 0; j < (n); j++)
//		{
//			h[j] = (double *)malloc((n) * sizeof(double));
//			for (int k = 0; k < (n); k++)
//			{
//				if (j == k)
//				{
//					h[j][j] = 1;
//				}
//				else
//				{
//					h[j][k] = 0;
//				}
//			}
//		}
//		for (int j = 1 + i; j < n; j++)
//		{
//			for (int k = 1 + i; k < n; k++)
//			{
//				h[j][k] -= p[j - i - 1][k - i - 1] * 2;
//			}
//		}
//		double **r = (double **)malloc((n) * sizeof(double *));
//		for (int j = 0; j < (n); j++)
//		{
//			r[j] = (double *)malloc((n) * sizeof(double));
//		}
//		for (int k = 0; k < n; k++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				double cur = 0;
//				for (int l = 0; l < n; l++)
//				{
//					double a = h[j][l] * A[l][k];
//					cur += a;
//				}
//				if (fabs(cur) < 1e-12)
//					cur = 0;
//				r[j][k] = cur;
//			}
//		}
//		for (int k = 0; k < n; k++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				A[k][j] = r[k][j];
//			}
//		}
//		for (int k = 0; k < n; k++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				double cur = 0;
//				for (int l = 0; l < n; l++)
//				{
//					cur += A[j][l] * h[l][k];
//				}
//				if (fabs(cur) < 1e-12)
//					cur = 0;
//				r[j][k] = cur;
//			}
//		}
//		for (int k = 0; k < n; k++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				A[k][j] = r[k][j];
//			}
//		}
//	}
//}
//
//void householder(double **A, double **Q, double **R, int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			R[i][j] = A[i][j];
//			if (i == j)
//			{
//				Q[i][j] = 1;
//			}
//			else
//			{
//				Q[i][j] = 0;
//			}
//		}
//	}
//
//	for (int k = 0; k < n - 1; k++)
//	{
//		double x_norm = 0;
//		for (int i = k; i < n; i++)
//		{
//			x_norm += R[i][k] * R[i][k];
//		}
//		x_norm = sqrt(x_norm);
//
//		double v[n];
//		for (int i = 0; i < n; i++)
//		{
//			if (i < k)
//			{
//				v[i] = 0;
//			}
//			else if (i == k)
//			{
//				v[i] = R[i][k] + copysign(x_norm, R[i][k]);
//			}
//			else
//			{
//				v[i] = R[i][k];
//			}
//		}
//
//		double v_norm = 0;
//		for (int i = k; i < n; i++)
//		{
//			v_norm += v[i] * v[i];
//		}
//		v_norm = sqrt(v_norm);
//
//		for (int i = k; i < n; i++)
//		{
//			v[i] /= v_norm;
//		}
//
//		for (int i = k; i < n; i++)
//		{
//			double temp = 0;
//			for (int j = k; j < n; j++)
//			{
//				temp += R[j][i] * v[j];
//			}
//			for (int j = k; j < n; j++)
//			{
//				R[j][i] -= 2 * temp * v[j];
//			}
//		}
//
//		for (int i = 0; i < n; i++)
//		{
//			double temp = 0;
//			for (int j = k; j < n; j++)
//			{
//				temp += Q[i][j] * v[j];
//			}
//			for (int j = k; j < n; j++)
//			{
//				Q[i][j] -= 2 * temp * v[j];
//			}
//		}
//	}
//}





// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
//
// void jacobi(int n, double A[], double eigenvalues[], double eigenvectors[], int *num_iterations);
//
// int main(int argc, char *argv[])
//{
//	int n;
//	FILE *f;
//	f = fopen(argv[1], "r");
//	fscanf(f, "%d", &n);
//	double *A = (double *)malloc(n * n * sizeof(double));
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			fscanf(f, "%lf", &A[i * n + j]);
//		}
//	}
//	fclose(f);
//
//	double *eigenvalues = (double *)malloc(n * sizeof(double));
//	double *eigenvectors = (double *)malloc(n * n * sizeof(double));
//
//	int num_iterations = 0;
//	jacobi(n, A, eigenvalues, eigenvectors, &num_iterations);
//
//	f = fopen(argv[2], "w");
//	for (int i = 0; i < n; i++)
//	{
//		fprintf(f, "%g\n", eigenvalues[i]);
//	}
//	fclose(f);
//	free(A);
//	free(eigenvalues);
//	free(eigenvectors);
//	return 0;
// }
//
// void jacobi(int n, double A[], double eigenvalues[], double eigenvectors[], int *num_iterations)
//{
//	double epsilon = 1e-10;	   // convergence criterion
//	double max_offdiag = 0.0;
//	int p, q;
//	double c, s;
//	double tau, t;
//
//	// Initialize eigenvectors to identity matrix
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			eigenvectors[i * n + j] = (i == j) ? 1.0 : 0.0;
//		}
//	}
//
//	// Iterate until convergence
//	while (1)
//	{
//		// Find maximum off-diagonal element
//		max_offdiag = 0.0;
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = i + 1; j < n; j++)
//			{
//				if (fabs(A[i * n + j]) > max_offdiag)
//				{
//					max_offdiag = fabs(A[i * n + j]);
//					p = i;
//					q = j;
//				}
//			}
//		}
//
//		if (max_offdiag < epsilon)
//		{
//			break;
//		}
//
//		// Compute Jacobi rotation matrix
//		tau = (A[q * n + q] - A[p * n + p]) / (2.0 * A[p * n + q]);
//		t = 1.0 / (fabs(tau) + sqrt(tau * tau + 1.0));
//		if (tau < 0.0)
//		{
//			t = -t;
//		}
//		c = 1.0 / sqrt(t * t + 1.0);
//		s = t * c;
//
//		// Update matrix A
//		double Apq = A[p * n + q];
//		A[p * n + q] = 0.0;
//		A[q * n + p] = 0.0;
//		A[p * n + p] -= t * Apq;
//		A[q * n + q] += t * Apq;
//		for (int k = 0; k < n; k++)
//		{
//			if (k != p && k != q)
//			{
//				double Akp = A[k * n + p];
//				double Akq = A[k * n + q];
//				A[k * n + p] = c * Akp - s * Akq;
//				A[p * n + k] = A[k * n + p];
//				A[k * n + q] = c * Akq + s * Akp;
//				A[q * n + k] = A[k * n + q];
//			}
//		}
//
//		// Update eigenvectors
//		for (int k = 0; k < n; k++)
//		{
//			double vkp = eigenvectors[k * n + p];
//			double vkq = eigenvectors[k * n + q];
//			eigenvectors[k * n + p] = c * vkp - s * vkq;
//			eigenvectors[k * n + q] = c * vkq + s * vkp;
//		}
//
//		(*num_iterations)++;
//	}
//
//	// Set eigenvalues
//	for (int i = 0; i < n; i++)
//	{
//		eigenvalues[i] = A[i * n + i];
//	}
// }
//
// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
//
// void jacobi(int n, double A[][n], double eigenvalues[], double eigenvectors[][n], int *num_iterations);
//
// int main(int argc, char *argv[])
//{
//	int n;
//	FILE *f;
//	f = fopen(argv[1], "r");
//	fscanf(f, "%d", &n);
//	double A[n][n];
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			fscanf(f, "%lf", &A[i][j]);
//		}
//	}
//	fclose(f);
//	double eigenvalues[n];
//	double eigenvectors[n][n];
//	int num_iterations = 0;
//	jacobi(n, A, eigenvalues, eigenvectors, &num_iterations);
//
//	f = fopen(argv[2], "w");
//	for (int i = 0; i < n; i++)
//	{
//		fprintf(f, "%g\n", eigenvalues[i]);
//	}
//	fclose(f);
//	return 0;
// }
//
// void jacobi(int n, double A[][n], double eigenvalues[], double eigenvectors[][n], int *num_iterations)
//{
//	double epsilon = 1e-10;	   // convergence criterion
//	double max_offdiag = 0.0;
//	int p, q;
//	double c, s;
//	double tau, t;
//
//	// Initialize eigenvectors to identity matrix
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			eigenvectors[i][j] = 0.0;
//		}
//		eigenvectors[i][i] = 1.0;
//	}
//
//	// Iterate until convergence
//	while (1)
//	{
//		// Find maximum off-diagonal element
//		max_offdiag = 0.0;
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = i + 1; j < n; j++)
//			{
//				if (fabs(A[i][j]) > max_offdiag)
//				{
//					max_offdiag = fabs(A[i][j]);
//					p = i;
//					q = j;
//				}
//			}
//		}
//
//		if (max_offdiag < epsilon)
//		{
//			break;
//		}
//
//		// Compute Jacobi rotation matrix
//		tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
//		t = 1.0 / (fabs(tau) + sqrt(tau * tau + 1.0));
//		if (tau < 0.0)
//		{
//			t = -t;
//		}
//		c = 1.0 / sqrt(t * t + 1.0);
//		s = t * c;
//
//		// Update matrix A
//		double Apq = A[p][q];
//		A[p][q] = 0.0;
//		A[q][p] = 0.0;
//		A[p][p] -= t * Apq;
//		A[q][q] += t * Apq;
//		for (int k = 0; k < n; k++)
//		{
//			if (k != p && k != q)
//			{
//				double Akp = A[k][p];
//				double Akq = A[k][q];
//				A[k][p] = c * Akp - s * Akq;
//				A[p][k] = A[k][p];
//				A[k][q] = c * Akq + s * Akp;
//				A[q][k] = A[k][q];
//			}
//		}
//
//		// Update eigenvectors
//		for (int k = 0; k < n; k++)
//		{
//			double vkp = eigenvectors[k][p];
//			double vkq = eigenvectors[k][q];
//			eigenvectors[k][p] = c * vkp - s * vkq;
//			eigenvectors[k][q] = c * vkq + s * vkp;
//		}
//
//		(*num_iterations)++;
//	}
//
//	// Set eigenvalues
//	for (int i = 0; i < n; i++)
//	{
//		eigenvalues[i] = A[i][i];
//	}
// }