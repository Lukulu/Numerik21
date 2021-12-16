/*@author: Yuliya Kryvistkaya, Jiani Shen, Laura Zywietz Rolón*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/*compile command:
gcc -Wall -pedantic -Wextra -Wno-unused-parameter H3.c -o ex -lm */

/*--------------------------------functions-----------------------------------*/
/*print the matrix A with
n: number of rows
m: number of columns*/
void printMatrix(int n, int m, double* A) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (j == 0) printf("( ");
			printf("%3.3lf ", A[j * n + i]);
			if (j == (m - 1)) printf(")\n");
		}
	}
	printf("\n");
}

/*print a vector v with length length
for debugging purpose*/
void printVector(int length, double* v) {
	// printf("The vector has the form:\n");
	for (int i = 0; i < length; i++) {
		if (i == 0) printf("( ");
		printf("%3.3lf ", v[i]);
		if (i == (length - 1)) printf(")\n");
	}
}

/*perform matrix multiplication of nxm matrix A and
mxl matrix B and store the result in nxl matrix AB*/
void matrix_matrix_mul(double* A, double* B, double* AB, int n, int m, int l){
	/*initialize result matrix AB to zero matrix*/
	for (int j = 0; j < l; j++)
	{
		for (int i = 0; i < n; i++)
		{
			AB[j*n+i] =0;
		}
	}

	/*perform matrix multiplication AB = A*B*/
	for (int j = 0; j < l; j++)
	{
		for (int i = 0; i < n; i++)
		{
			for (int d = 0; d < m; d++)
			{
				AB[j*n+i] += A[d*m+i]*B[j*l+d];
			}
		}
	}
}

/*perform matrix vector multiplication of nxn matrix A and
vector x of length n and store the result in the vector Ax*/
void matrix_vector_mul(double* A, double* x, double* Ax, int n){
	/*initialize result vector Ax to zero vector*/
	for (int i = 0; i < n; i++)
	{
		Ax[i] = 0;
	}

	/*perform matrix-vector multiplication Ax = A*x*/
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Ax[i] += A[j*n+i]*x[j];
		}
	}
}

/*caculate euclidean scalar product of two vectors vec1 and vec2
 of length n*/
double euclidean_scp(double* vec1, double* vec2, int length){
	double res = 0;
	for (int i = 0; i < length; i++) {
		res += vec1[i] * vec2[i];
	}
	return res;
}

/*print the transposed mxn matrix At of nxm matrix A*/
void transpose_matrix(double *A, double* At, int n, int m){
	/*initialize result matrix AB to zero matrix*/
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
		{
			At[j*n+i] =0;
		}
	}

 for (int i = 0; i < m; i++)
 {
	 for (int j = 0; j < n; j++)
	 {
		 At[j*n+i] = A[i*n+j];
	 }
 }
}

/*calculate (R,R) and (P, AP) and print it to the command line*/
void print_results_algorithms(double* R, double* P, double* A, int n, int K){
	double* Rt = (double*)malloc(sizeof(double) * K * n);
	transpose_matrix(R, Rt, n, K);

	double* Pt = (double*)malloc(sizeof(double) * K * n);
	transpose_matrix(P, Pt, n, K);

	double* AP = (double*)malloc(sizeof(double) * n * K);
	double* RtR = (double*)malloc(sizeof(double) * K * K);
	double* PtAP = (double*)malloc(sizeof(double) * K * K);

	matrix_matrix_mul(Rt, R, RtR, K, n, K);
	matrix_matrix_mul(A, P, AP, n, n, K);
	matrix_matrix_mul(Pt, AP, PtAP, K, n, K);

	printf("The scalar product of matrices (R,R) has the form: \n");
	printMatrix(K, K, RtR);
	printf("The scalar product of matrices (P,AP) has the form: \n");
	printMatrix(K, K, PtAP);

	free(Rt);
	free(Pt);
	free(AP);
	free(RtR);
	free(PtAP);
}
/*-------------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
	double* A; //input matrix A
	int n = 0; //matrix and vector dimension
	int K = 10; //maximal number of iteration steps
	double* b; //input vector b
	double* x; //solution of the LSE
	double eta = 1e-8; //precision of approximation

	/*----------------------Read Input LSE from Command Line-----------------------*/
	printf("We want to solve the system of linear equations of the form Ax=b:\n");

	/*Read number of rows n, number of columns m and the matrix A itself
	from console and store them in the corresponding containers.*/
	printf("Please enter the dimension of the system: ");
	scanf("%d", &n);

	//terminate program if n = 0
	if (n == 0)
	{
		printf("Error! The number of rows and columns has to be different from zero.\n");
		return 1;
	}

	/*the matrix is columnwise read in to allow for a
	faster memory access of its columns.*/
	printf("Please enter the lower diagonal of the matrix A columnwise: ");
	A = (double*)malloc(sizeof(double) * n * n);
	for (int j = 0; j < n; j++) //iterate over columns
	{
		for (int i = j; i < n; i++) //iterate over rows
		{
			scanf("%lf", &A[j * n + i]);
      if(i!=j) A[i*n+ j] = A[j * n + i];
		}
	}

	/*Read vector b from the command line*/
	printf("Please enter your vector b: ");
	b = (double*)malloc(sizeof(double) * n);
	for (int i = 0; i < n; i++)
	{
		scanf("%lf", &b[i]);
	}

	/*Read vector x from the command line*/
	printf("Please enter starting value of your solution vector x: ");
	x = (double*)malloc(sizeof(double) * n);
	for (int i = 0; i < n; i++)
	{
		scanf("%lf", &x[i]);
	}

	/*Print inputs saved from the command line*/
	printf("Your input matrix A has the form: \n");
	printMatrix(n, n, A);

	printf("Your input vector b has transposed the form: \n");
	printVector(n, b);

	printf("\nYour input vector x has transposed the form: \n");
	printVector(n, x);

	/*select algorithm from input in the commeand line*/
	int alg_num = 0;
	printf("\nPlease choose the algorithm: \n");
	printf("Steepest Decent (SD(A,b,x,eta)): 1\n");
	printf("Conjugated Gradient (CG(A,b,x,eta)): 2\n");
	scanf("%d", &alg_num);
	if(alg_num != 1 && alg_num!=2){
		printf("Error, no valid number was selected. Exit\n");
		return 1;
	}

	/*---------------------------algorithms-----------------------------------*/
	/*define matrix P of search directions: P = (p_1,p_2,...,p_K)
	and the matrix R of residues: R = (r_1,r_2,...,r_K)*/
	double *P;
	double *R;
	P = (double*)malloc(sizeof(double) * n * K);
	R = (double*)malloc(sizeof(double) * n * K);
	for (int j = 0; j < K; j++)
	{
		for (int i = 0; i < n; i++)
		{
			P[j*n+i] = 0;
			R[j*n+i] = 0;
		}
	}

	/*----------------------Steepest Decent ----------------------------------*/
	if(alg_num == 1)
	{
		printf("Executing Steepest Decent...\n");
		double norm_r = 0, norm_b = 0;

		double* Ax = (double*)malloc(sizeof(double)*n);
		matrix_vector_mul(A, x, Ax, n);

		double* r = (double*)malloc(sizeof(double)*n);
		for (int i = 0; i < n; i++) r[i] = b[i] - Ax[i];
		free(Ax);

		/*calculate norm of r*/
		norm_r = euclidean_scp(r, r, n);
		norm_r = sqrt(norm_r);

		/*calculate norm of b*/
		norm_b = euclidean_scp(b, b, n);
		norm_b = sqrt(norm_b);

		int k = 0;
		while (norm_r > eta*norm_b && k < K)
		{
			// printf("\nIterations step: %d\n", k);
			// printf("------------------------\n");
			// printf("eta*norm: %e\n", eta*norm_b);
			// printf("norm_r: %e\n", norm_r);
			double* p = (double*)malloc(sizeof(double)*n);
			matrix_vector_mul(A, r, p, n);

			double alpha = euclidean_scp(r,r,n) / euclidean_scp(p,r,n);

			for (int i = 0; i < n; i++)
			{
				x[i] += alpha*r[i];
				r[i] -= alpha*p[i];
				P[k*n+i] = p[i];
				R[k*n+i] = r[i];
			}
			//adjust norm
			norm_r = euclidean_scp(r, r, n);
			norm_r = sqrt(norm_r);

			k++;
			free(p);
		}
		free(r);

		/*Output*/
		printf("Solution for x:\n");
		printVector(n, x);

		printf("Number of iteration steps: %d\n", k);

		/*calculate (R,R) and (P, AP) and print it to the command line*/
		// print_results_algorithms(R, P, A, n, K);
	}
	/*----------------------Conjugate Gradient--------------------------------*/
  // /*set x to zero*/
	// for (int i = 0; i < n; i++) {
	// 	x[i] = 0;
	// }

	if(alg_num == 2)
	{
		printf("Executing CG...\n");
		double norm_r = 0, norm_b = 0;

		/*initialize r (residual) and x for the CG algorithm with r = b and x = 0*/
		double* new_r = (double*)malloc(sizeof(double)*n);
		double* r = (double*)malloc(sizeof(double)*n);
		double* p = (double*)malloc(sizeof(double)*n);
		double* Ax = (double*)malloc(sizeof(double)*n);
		matrix_vector_mul(A, x, Ax, n);

		for (int i = 0; i < n; i++) {
			r[i] = b[i]- Ax[i];
			p[i] = r[i];
			new_r[i] = r[i];
		}
		free(Ax);

		/*calculate norm of r*/
		norm_r = euclidean_scp(r, r, n);
		norm_r = sqrt(norm_r);

		/*calculate norm of b*/
		norm_b = euclidean_scp(b, b, n);
		norm_b = sqrt(norm_b);

		int k = 0;
		while (norm_r > eta*norm_b && k < K)
		{
			// printf("\nIterations step: %d\n", k);
			// printf("------------------------\n");
			// printf("eta*norm: %e\n", eta*norm_b);
			// printf("norm_r: %e\n", norm_r);
			double* Ap = (double*)malloc(sizeof(double)*n);
			matrix_vector_mul(A, p, Ap, n);

			double alpha = euclidean_scp(r,r,n) / euclidean_scp(p,Ap,n);

			for (int i = 0; i < n; i++)
			{
				x[i] += alpha*p[i];
				new_r[i] = r[i] - alpha*Ap[i];
			}
			double beta = euclidean_scp(new_r,new_r,n)/euclidean_scp(r,r,n);

			for (int i = 0; i < n; i++)
			{
				p[i] = new_r[i] + beta*p[i];
				P[k*n+i] = p[i];
				R[k*n+i] = new_r[i];
				r[i] = new_r[i];
			}
			/*adjust norm of r*/
			norm_r = euclidean_scp(new_r, new_r, n);
			norm_r = sqrt(norm_r);

			k++;
			free(Ap);
		}
		free(new_r);
		free(r);
		free(p);

		printf("Solution for x:\n");
		printVector(n, x);

		printf("Number of iteration steps: %d\n", k);

		/*calculate (R,R) and (P, AP) and print it to the command line*/
		// print_results_algorithms(R, P, A, n, K);
	}
	/*--------------------------------------------------------------------------*/

	/*free allocated memory*/
  free(A);
	free(R);
	free(P);
	free(x);
	free(b);

	return 0;
}
