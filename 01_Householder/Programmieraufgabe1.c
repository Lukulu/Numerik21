//@author: Laura Zywietz Rolón, Alina Becker
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/*compile command:
gcc -Wall -pedantic -Wextra -Wno-unused-parameter H1.c -o h1 -lm */

/*--------------------------------functions-----------------------------------*/
/*print the matrix A with
n: number of rows
m: number of columns*/
void printMatrix(int n, int m, double *A){
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			printf("%3.3lf ", A[j*n+i]);
			if (j==(m-1)){ printf("\n");}
		}
	}
	printf("\n");
}

 /*print the transposed matrix A*/
void printTransposedMatrix(int n, int m, double *A){
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			printf("%3.3lf ", A[i*n+j]);
			if (j==(m-1)){ printf("\n");}
		}
	}
	printf("\n");
}

/*print a vector v with length length
for debugging purpose*/
void printVector(int length, double*v){
	printf("The vector has the form:\n");
	for (int i = 0; i < length; i++) {
		if(i==0) printf("( ");
		printf("%3.3lf ", v[i]);
		if(i==(length-1)) printf(")\n");
	}
}

/*calculate minimum of two numbers a, b*/
int findMin(int a, int b){
	if (a < b) return a;
	else return b;
}

/*determine the sign of a given number a*/
int sgn(double a){
	if(a >= 0) return 1;
	else return -1;
}

/*determine absolute value of column vector at pos of A
of length n according to the Euclidean norm.
If squared = 1: return the squared value to avoid
its calculation afterwars by squaring the root.*/
double abs_vec(double *A, int pos, int n, int squared){
	double res = 0;
	for (int i = pos; i < n; i++) {
		res += A[pos*n+i]*A[pos*n+i];
	}
	if(squared) return res;
	else return sqrt(res);
}

/*determine the absolute value of a vector v of length length-pos.
If squared = 1: return the squared value to avoid
its calculation afterwars by squaring the root.*/
double abs_vec1(double *v, int pos, int length, int squared){
	double res = 0;
	for (int i = pos; i < length; i++) {
		res += v[i]*v[i];
	}
	if(squared) return res;
	else return sqrt(res);
}

/*calculate scalar product of two vectors a and b of length length*/
double skp(double *a, double *b, int length){
	double res = 0;
	for (int i = 0; i < length; i++) {
		res += a[i]*b[i];
	}
	return res;
}

// /*calculate absolute value*/
// double abs_vec2(double *A, int pos, int n){
// 	double res = 0;
// 	for (int i = pos; i < n; i++) {
// 		res += A[pos*n+i]*A[pos*n+i];
// 	}
// 	return sqrt(res);
// }

/*calculate product of matrix Q with given vector w of dimension n,
where v depicts the HH-vector of iteration step j,
A the input matrix A */
void Q_times_w(double *w, double *v, double *A, int n, int j)
{
	double *tmp_v = (double *) malloc(sizeof(double)*n);
	/*expand the HH-vector to a vector of dimension n*/
	for (int i = 0; i < n; i++)
	{
		if (i < j) tmp_v[i] = 0;
		else tmp_v[i] = v[i]/v[j];
	}
	printVector(n, tmp_v);

	double factor = 1. + fabs(A[j*n+j])/abs_vec1(A, j, n, 0);
	double skp_factor = skp(tmp_v, w, n);
	//overwrite vector w with result of computation
	for (int i = 0; i < n; i++)
	{
		w[i] = w[i] - factor*skp_factor*tmp_v[i];
	}
	free(tmp_v);
}
/*----------------------------------------------------------------------------*/

int main (int argc, char** argv)
{
	double *A; //input matrix A
	double *v; //HH vector of current iteration step j
	double *x; //j-th column of A
	int n = 0, m = 0; //number of rows/columns

	double *Q; //matrix Q
	double *w; //store 1,...,m unit vectors of R^n temporarily

	/*----------------------Read Matrix from Command Line-----------------------*/
	/*Read number of rows n, number of columns m and the matrix A itself
	from console and store them in the corresponding containers.*/
	printf("Please enter the matrix size:\n");
	printf("Number of rows: ");
	scanf("%d", &n);
	printf("Number of columns: ");
	scanf("%d", &m);

	//terminate program if n = 0 or m = 0
	if(n == 0 || m == 0)
	{
		printf("Error! The number of rows and columns has to be different from zero. Try again.\n");
		return 1;
	}

	/*the matrix is columnwise read in to allow for a
	faster memory access of its columns.*/
	printf("Please enter the matrix column-wise:");
	A = (double *) malloc(sizeof(double) * n* m);
	for (int j = 0; j < m; j++) //iterate over columns
	{
		for (int i = 0; i < n; i++) //iterate over rows
		{
			scanf("%lf", &A[j*n+i]);
		}
	}

	// n = 3;
	// m = 3;
	// A = (double *) malloc(sizeof(double) * n* m);
	// // A[0] = 2;
	// // A[1] = 1;
	// // A[2] = 2;
	// // A[3] = 4;
	// // A[4] = 1;
	// // A[5] = -3;
	// // A[6] = -4;
	// // A[7] = 2;
	// // A[8] = 0;
	// A[0] = 1;
	// A[1] = 2;
	// A[2] = 3;
	// A[3] = 4;
	// A[4] = 5;
	// A[5] = 6;
	// A[6] = 7;
	// A[7] = 8;
	// A[8] = 9;

	printf("Your input matrix has the form: \n");
	printMatrix(n, m, A);
	/*--------------------------------------------------------------------------*/

	/*---------------------Calculate QR Decomposition of A----------------------*/
	/*allocate nxm-matrix Q and fill it with unit vectors of R^m*/
	Q = (double *) malloc(sizeof(double)*n*m);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			Q[j*n+i] = 0;
			if(j == i) Q[j*n+i] = 1;
		}
	}
	/*allocate vector w of size n, and set it to the first column of Q*/
	w = (double *) malloc(sizeof(double)*n);
	for (int i = 0; i < n; i++) w[i] = Q[i];

	v = (double *) malloc(sizeof(double)*n); //buffer Householder vector of current iteration step
	x = (double *) malloc(sizeof(double)*n); //buffer column vector of A at current iteration step
	//set x to the first column of A to start the iteration
	for (int i = 0; i < n; i++) x[i] = A[i];

	/*perform QR-decomposition: calculate HH-vectors and R*/
	for (int j = 0; j <= findMin(m-1, n); j++) //min(m-1, n)+1 iteration steps
	{
		printf("Iterationstep: %d\n", j);

		double abs = abs_vec(A,j,n, 0); //absolute value of column vector j of A
		// printf("Absoulte value of column vector j of A: %2.1lf \n", abs);

		/*calculate Householder vector v for iteration step j*/
		for (int i = j; i < n; i++)
		{
			if(i==j) v[i] = A[j*n+i] + sgn(A[j*n+i])*abs;
			else     v[i] = A[j*n+i];
		}
		printVector(n, v);

		/*calculate Q*/
		for (int l = 0; l < m; l++) //iterate over all column vectors of Q
		{
			Q_times_w(w, v, A, n, j); //calculate l-th column of Q
			for (int k = j; k < n; k++) Q[l*n+k] = w[k]; //store result in Q
			for (int k = 0; k < n; k++) w[k] = Q[(l+1)*n+k]; //reset w to next column of Q
		}
		for (int k = 0; k < n; k++) w[k] = Q[k]; //reset w to the first column of Q

		// //calculate Q
		// for(int l = 0; l < m; l++) //iterate over columns
		// {
		// 	for(int i = j; i < n; i++) w[i] = Q[i*n+l]; //transposed vector
		//
		// 	for(int k = j; k < n; k++)
		// 	{
		// 		for (int d = j; d < n; d++) Q[k*n+l] += -2/abs_v*v[k]*v[d]*w[d];
		// 	}
		// }

		//calculate H_j * A_j
		double abs_v = abs_vec1(v,j,n, 1); //absolute value of HH vector v squared
		// printf("Absolute value HH vector v: %2.1lf\n", abs_v);
		for (int l = j; l < m; l++) //iterate over columns
		{
			for (int k = j; k < n; k++) //iterate over rows
			{
				for (int d = j; d < n; d++)
				{
					A[l*n+k] += -2/abs_v*v[k]*v[d]*x[d];
				}
			}
			/*set vector x to the next column of A*/
			for (int k = 0; k < n; k++) x[k] = A[(l+1)*n+k];
		}

		/*reset vector x to the first column of
		the new matrix A*/
		for (int k = 0; k < n; k++) x[k] = A[(j+1)*n+k];

		/*write the normalized HH vector into the
		j-th column of A below the diagonal.*/
		if(j!=findMin(m-1, n))
		{
			for (int k = (j+1); k < n; k++) A[j*n+k] = v[k]/v[j];
		}

		printf("\nResulting matrix of the iteration step %d with HH vectors: \n", j);
		printMatrix(n, m, A);

		printf("\nResulting matrix Q of the iteration step %d: \n", j);
		printTransposedMatrix(n, m, Q);
	}
	/*--------------------------------------------------------------------------*/

	/*free allocated memory*/
	free(A);
	free(v);
	free(x);
	free(Q);
	free(w);

	return 0;
}
