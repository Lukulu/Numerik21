/**/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <stdbool.h>

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
	printf("The vector has the form:\n");
	for (int i = 0; i < length; i++) {
		if (i == 0) printf("( ");
		printf("%3.3lf ", v[i]);
		if (i == (length - 1)) printf(")\n");
	}
}

/*perform matrix multiplication of nxm matrix A and
mxn matrix B and store the result in nxn matrix AB*/
void matrix_matrix_mul(double* A, double* B, double* AB, int n, int m){
	/*initialize result matrix AB to zero matrix*/
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			AB[j*n+i] =0;
		}
	}

	/*perform matrix multiplication AB = A*B*/
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			for (int d = 0; d < m; d++)
			{
				AB[j*n+i] += A[d*m+i]*B[j*n+d];
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
		Ax[i] =0;
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

double euclidean_scp(double* vec1, double* vec2, int length){
	double res = 0;
	for (int i = 0; i < length; i++) {
		res += vec1[i] * vec2[i];
	}
	return res;
}

// /*determine the sign of a given number a*/
// int sgn(double a) {
// 	if (a >= 0) return 1;
// 	else return -1;
// }
//
// /*calculate minimum of two numbers a, b*/
// int findMin(int a, int b) {
// 	if (a < b) return a;
// 	else return b;
// }

/*determine absolute value of column vector at pos of A
of length n according to the Euclidean norm.
If squared = 1: return the squared value to avoid
its calculation afterwars by squaring the root.*/
double abs_vec(double* A, int pos, int n, int squared) {
	double res = 0;
	for (int i = pos; i < n; i++) {
		res += A[pos * n + i] * A[pos * n + i];
	}
	if (squared) return res;
	else return sqrt(res);
}

double abs_vec2(double* A, int pos, int n, int squared) {
	double res = 0;
	for (int i = pos; i < n-1; i++) {
		res += A[pos * n + i+1] * A[pos * n + i+1];
		// printf("i=%d: %lf\n", i, A[pos * n + i+1]);
	}
	// printf("\n");
	if (squared) return res;
	else return sqrt(res);
}

/*determine the absolute value of a vector v of length length-pos.
If squared = 1: return the squared value to avoid
its calculation afterwars by squaring the root.*/
double abs_vec1(double* v, int pos, int length, int squared) {
	double res = 0;
	for (int i = pos; i < length; i++) {
		res += v[i] * v[i];
	}
	if (squared) return res;
	else return sqrt(res);
}

// /*calculate QR-decomposition of a matrix A of dimension nxm*/
// void get_QR_decomposition(double* A, double* Q, int n, int m){
//   /*allocate vector w of size n, and set it to the first column of Q*/
//   double* w = (double*)malloc(sizeof(double) * n);
//   for (int i = 0; i < n; i++) w[i] = Q[i];
//
//   /*buffer column vector of A at current iteration step*/
//   double* x = (double*)malloc(sizeof(double) * n);
//   /*set x to the first column of A to start the iteration*/
//   for (int i = 0; i < n; i++) x[i] = A[i];
//
//   /*buffer Householder vector of current iteration step*/
//   double* v = (double*)malloc(sizeof(double) * n);
//
//   /*perform QR-decomposition: calculate HH-vectors and R*/
//   for (int j = 0; j <= findMin(m - 1, n); j++) //min(m-1, n)+1 iteration steps
//   {
//     double abs = abs_vec(A, j, n, 0); //absolute value of column vector j of A
//
//     /*calculate Householder vector v for iteration step j*/
//     for (int i = j; i < n; i++)
//     {
//       if (i == j) v[i] = A[j * n + i] + sgn(A[j * n + i]) * abs;
//       else     v[i] = A[j * n + i];
//     }
//
//     //calculate H_j * A_j
//     double abs_v = abs_vec1(v, j, n, 1); //absolute value of HH vector v squared
//     for (int l = j; l < m; l++) //iterate over columns
//     {
//       for (int k = j; k < n; k++) //iterate over rows
//       {
//         for (int d = j; d < n; d++)
//         {
//           A[l * n + k] += -2 / abs_v * v[k] * v[d] * x[d];
//         }
//       }
//       /*set vector x to the next column of A*/
//       for (int k = 0; k < n; k++) x[k] = A[(l + 1) * n + k];
//     }
//
//     /*calculate matrix Q for each iteration step j*/
//     for (int l = 0; l < n; l++) //iterate over columns of Q
//     {
//       //set vector w to transposed column vector
//       for (int i = j; i < n; i++) w[i] = Q[i * n + l];
//       /*calculate next Q-matrix analogeously to part 1 above, taking the transposition
//       for the final result already into account*/
//       for (int k = j; k < n; k++)
//       {
//         for (int d = j; d < n; d++) Q[k * n + l] += -2 / abs_v * v[k] * v[d] * w[d];
//       }
//     }
//
//     /*reset vector x to the first column of the new matrix A*/
//     for (int k = 0; k < n; k++) x[k] = A[(j + 1) * n + k];
//
//     // /*write the normalized HH vector into thej-th column of A below the diagonal.*/
//     // if (j != findMin(m - 1, n))
//     // {
//     //   for (int k = (j + 1); k < n; k++) A[j * n + k] = v[k] / v[j];
//     // }
//   }
//   free(w);
//   free(x);
//   free(v);
// }
//
// /*calculate QR-decomposition of a band matrix A of dimension nxn with*/
// void get_QR_decomp_band_matrix(double* A, double* Q, double* R, int n){
//   /*allocate vector w of size n, and set it to the first column of Q*/
//   double* w = (double*)malloc(sizeof(double) * n);
//   for (int i = 0; i < n; i++) w[i] = Q[i];
//
//   /*buffer column vector of A at current iteration step*/
//   double* x = (double*)malloc(sizeof(double) * n);
//   /*set x to the first column of A to start the iteration*/
//   for (int i = 0; i < n; i++) x[i] = A[i];
//
//   /*buffer Householder vector of current iteration step*/
//   double* v = (double*)malloc(sizeof(double) * n);
//
// 	/*set R = A and Q = I to start the iteration*/
// 	for (int j = 0; j < n; j++) {
// 		for (int i = 0; i < n; i++) {
// 			R[j*n+i] = A[j*n+i];
// 			Q[j * n + i] = 0;
// 			if (j == i) Q[j * n + i] = 1;
// 		}
// 	}
//
//   /*perform QR-decomposition: calculate HH-vectors and R*/
//   for (int j = 0; j < n; j++) //n iteration steps
//   {
//     double abs = abs_vec(R, j, n, 0); //absolute value of column vector j of A
//
//     /*calculate Householder vector v for iteration step j*/
//     for (int i = j; i < n; i++)
//     {
//       if (i == j) v[i] = R[j * n + i] + sgn(R[j * n + i]) * abs;
//       else     v[i] = R[j * n + i];
//     }
//
//     //calculate H_j * A_j
//     double abs_v = abs_vec1(v, j, n, 1); //absolute value of HH vector v squared
//     for (int l = j; l < n; l++) //iterate over columns
//     {
//       for (int k = j; k < n; k++) //iterate over rows
//       {
//         for (int d = j; d < n; d++)
//         {
//           R[l * n + k] += -2 / abs_v * v[k] * v[d] * x[d];
//         }
//       }
//       /*set vector x to the next column of A*/
//       for (int k = 0; k < n; k++) x[k] = R[(l + 1) * n + k];
//     }
//
//     /*calculate matrix Q for each iteration step j*/
//     for (int l = 0; l < n; l++) //iterate over columns of Q
//     {
//       //set vector w to transposed column vector
//       for (int i = j; i < n; i++) w[i] = Q[i * n + l];
//       /*calculate next Q-matrix analogeously to part 1 above, taking the transposition
//       for the final result already into account*/
//       for (int k = j; k < n; k++)
//       {
//         for (int d = j; d < n; d++) Q[k * n + l] += -2 / abs_v * v[k] * v[d] * w[d];
//       }
//     }
//     /*reset vector x to the first column of the new matrix A*/
//     for (int k = 0; k < n; k++) x[k] = R[(j + 1) * n + k];
//   }
//   free(w);
//   free(x);
//   free(v);
// }

/*print the transposed matrix At of A*/
void transposeMatrix(int n, double *A, double* At){
	/*initialize result matrix AB to zero matrix*/
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			At[j*n+i] =0;
		}
	}

 for (int i = 0; i < n; i++)
 {
	 for (int j = 0; j < n; j++)
	 {
		 At[j*n+i] = A[i*n+j];
	 }
 }
}
/*----------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
	double* A; //input matrix A
	int n = 0; //matrix and vector dimension
	int K = 1e5; //maximal number of iteration steps
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

	/*Read vector X from the command line*/
	printf("Please enter starting value of your solution vector x: ");
	x = (double*)malloc(sizeof(double) * n);
	for (int i = 0; i < n; i++)
	{
		scanf("%lf", &x[i]);
	}

	/*Print inputs saved from the command line*/
	printf("\nYour input matrix A has the form: \n");
	printMatrix(n, n, A);

	printf("Your input vector b has the form: \n");
	printVector(n, b);

	printf("Your input vector x has the form: \n");
	printVector(n, x);

	/*select algorithm from input in the commeand line*/
	int alg_num = 0;
	printf("Please choose the algorithm: \n");
	printf("Steepest Decent (SD(A,b,x,eta)): 1\n");
	printf("Conjugated Gradient (CG(A,b,x,eta)): 2\n");
	scanf("%d", &alg_num);
	if(alg_num != 1 && alg_num!=2){
		printf("Error, no valid number was selected. Exit\n");
		return 1;
	}

	/*test matrices*/
	// //A = Id(nxn)
	// for (int j = 0; j < n; j++) {
	// 	for (int i = 0; i < n; i++) {
	// 		A[j*n+i] = 0;
	// 		if (i==j) A[j*n+i] = 1;
	// 	}
	// }
	// //A = D(2^0, 2^1, ..., 2^(n-1))
	// for (int j = 0; j < n; j++) {
	// 	for (int i = 0; i < n; i++) {
	// 		A[j*n+i] = 0;
	// 		if (i==j) A[j*n+i] = pow(2, i);
	// 	}
	// }
	// //A_ij = 1/(i+j-1) (Hilbertmatrix)
	// for (int j = 0; j < n; j++) {
	// 	for (int i = 0; i < n; i++) {
	// 		A[j*n+i] = 1/(i+j+1);
	// 	}
	// }

	/*----------------------Steepest Decent ----------------------------------*/
	/*define matrix P of search directions: P = (p_1,p_2,...,p_K)
	and the matrix R of residues: R = (r_1,r_2,...,r_K)*/
	double *P;
	double *R;
	P = (double*)malloc(sizeof(double) * n * K);
	R = (double*)malloc(sizeof(double) * n * K);
	for (int j = 0; j < K; j++) {
		for (int i = 0; i < n; i++) {
			P[j*K+i] = 0;
			R[j*K+i] = 0;
		}
	}

	double* r;
	double* Ax;
	double* p;
	double norm_r = 0, norm_b = 0;

	Ax = (double*)malloc(sizeof(double)*n);
	matrix_vector_mul(A, x, Ax, n);

	r = (double*)malloc(sizeof(double)*n);
	for (int i = 0; i < n; i++) {
		r[i] = b[i] - Ax[i];
	}
	free(Ax);

	/*calculate norm of r*/
	norm_r = euclidean_scp(r, r, n);
	norm_r = sqrt(norm_r);

	/*calculate norm of b*/
	norm_b = euclidean_scp(b, b, n);
	norm_b = sqrt(norm_b);

	if(alg_num == 1)
	{
		int k = 0;
		while (norm_r > eta*norm_b || k <= K)
		{
			p = (double*)malloc(sizeof(double)*n);
			matrix_vector_mul(A, r, p, n);

			double alpha = euclidean_scp(r,r,n) / euclidean_scp(p,r,n);
			for (int i = 0; i < n; i++)
			{
				x[i] += alpha*r[i];
				r[i] -= alpha*p[i];
				P[k*K+i] = p[i];
				R[k*K+i] = r[i];
			}
			k++;
		}

		/*calculate (R,R) and (P, AP) and print it to the command line*/
		double* Rt = (double*)malloc(sizeof(double) * K * n);
		transposeMatrix(n, R, Rt);
		double* Pt = (double*)malloc(sizeof(double) * K * n);
		transposeMatrix(n, P, Pt);
		double* AP = (double*)malloc(sizeof(double) * n * K);
		double* RtR = (double*)malloc(sizeof(double) * K * K);
		double* PtAP = (double*)malloc(sizeof(double) * K * K);
		// matrix_matrix_mul(Rt, R, RtR, K, n);
		// matrix_matrix_mul(A, P, AP, n, K);
		// matrix_matrix_mul(Pt, AP, PtAP, n, n);

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
	/*----------------------Conjugate Gradient----------------------------------*/
	if(alg_num == 2)
	{
		int k = 0;
		while (norm_r > eta*norm_b)
		{
			p = (double*)malloc(sizeof(double)*n);
			matrix_vector_mul(A, r, p, n);

			double alpha = euclidean_scp(r,r,n) / euclidean_scp(p,r,n);
			for (int i = 0; i < n; i++)
			{
				x[i] += alpha*r[i];
				r[i] -= alpha*p[i];
				P[k*K+i] = p[i];
				R[k*K+i] = r[i];
			}
			k++;
		}

		/*calculate (R,R) and (P, AP) and print it to the command line*/
		double* Rt = (double*)malloc(sizeof(double) * K * n);
		transposeMatrix(n, R, Rt);
		double* Pt = (double*)malloc(sizeof(double) * K * n);
		transposeMatrix(n, P, Pt);
		double* AP = (double*)malloc(sizeof(double) * n * K);
		double* RtR = (double*)malloc(sizeof(double) * K * K);
		double* PtAP = (double*)malloc(sizeof(double) * K * K);
		// matrix_matrix_mul(Rt, R, RtR, K, n);
		// matrix_matrix_mul(A, P, AP, n, K);
		// matrix_matrix_mul(Pt, AP, PtAP, n, n);

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
	/*--------------------------------------------------------------------------*/

	/*free allocated memory*/
  free(A);
	free(R);
	free(P);
	free(x);
	free(b);
	free(p);
	free(r);

	return 0;
}
