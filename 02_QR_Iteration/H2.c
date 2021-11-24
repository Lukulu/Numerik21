/*@author: Laura Zywietz Rolon, Alina Becker*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/*compile command:
gcc -Wall -pedantic -Wextra -Wno-unused-parameter numerik21_exercise2.c -o ex2 -lm */

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

void matrix_mul(double* A, double* B, double* AB, int n){
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
			for (int d = 0; d < n; d++)
			{
				AB[j*n+i] += A[d*n+i]*B[j*n+d];
			}
		}
	}
}

/*calculate minimum of two numbers a, b*/
int findMin(int a, int b) {
	if (a < b) return a;
	else return b;
}

/*determine the sign of a given number a*/
int sgn(double a) {
	if (a >= 0) return 1;
	else return -1;
}

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

/*calculate QR-decomposition of a matrix A of dimension nxm*/
void get_QR_decomposition(double* A, double* Q, int n, int m){
  /*allocate vector w of size n, and set it to the first column of Q*/
  double* w = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; i++) w[i] = Q[i];

  /*buffer column vector of A at current iteration step*/
  double* x = (double*)malloc(sizeof(double) * n);
  /*set x to the first column of A to start the iteration*/
  for (int i = 0; i < n; i++) x[i] = A[i];

  /*buffer Householder vector of current iteration step*/
  double* v = (double*)malloc(sizeof(double) * n);

  /*perform QR-decomposition: calculate HH-vectors and R*/
  for (int j = 0; j <= findMin(m - 1, n); j++) //min(m-1, n)+1 iteration steps
  {
    double abs = abs_vec(A, j, n, 0); //absolute value of column vector j of A

    /*calculate Householder vector v for iteration step j*/
    for (int i = j; i < n; i++)
    {
      if (i == j) v[i] = A[j * n + i] + sgn(A[j * n + i]) * abs;
      else     v[i] = A[j * n + i];
    }

    //calculate H_j * A_j
    double abs_v = abs_vec1(v, j, n, 1); //absolute value of HH vector v squared
    for (int l = j; l < m; l++) //iterate over columns
    {
      for (int k = j; k < n; k++) //iterate over rows
      {
        for (int d = j; d < n; d++)
        {
          A[l * n + k] += -2 / abs_v * v[k] * v[d] * x[d];
        }
      }
      /*set vector x to the next column of A*/
      for (int k = 0; k < n; k++) x[k] = A[(l + 1) * n + k];
    }

    /*calculate matrix Q for each iteration step j*/
    for (int l = 0; l < n; l++) //iterate over columns of Q
    {
      //set vector w to transposed column vector
      for (int i = j; i < n; i++) w[i] = Q[i * n + l];
      /*calculate next Q-matrix analogeously to part 1 above, taking the transposition
      for the final result already into account*/
      for (int k = j; k < n; k++)
      {
        for (int d = j; d < n; d++) Q[k * n + l] += -2 / abs_v * v[k] * v[d] * w[d];
      }
    }

    /*reset vector x to the first column of the new matrix A*/
    for (int k = 0; k < n; k++) x[k] = A[(j + 1) * n + k];

    /*write the normalized HH vector into thej-th column of A below the diagonal.*/
    if (j != findMin(m - 1, n))
    {
      for (int k = (j + 1); k < n; k++) A[j * n + k] = v[k] / v[j];
    }
  }
  free(w);
  free(x);
  free(v);
}
/*----------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
	double* A; //input matrix A
	double* Q; //matrix Q
	int n = 0; //number of rows/columns
  // double eps = 0.000001 //check whether the HH-vector is zero, i.e. smaller than eps

	/*----------------------Read Matrix from Command Line-----------------------*/
	// /*Read number of rows n, number of columns m and the matrix A itself
	// from console and store them in the corresponding containers.*/
	// printf("Please enter the matrix size:\n");
	// printf("Number of rows and columns: ");
	// scanf("%d", &n);
	//
	// //terminate program if n = 0
	// if (n == 0)
	// {
	// 	printf("Error! The number of rows and columns has to be different from zero. Try again.\n");
	// 	return 1;
	// }
	//
	// /*the matrix is columnwise read in to allow for a
	// faster memory access of its columns.*/
	// printf("Please enter the lower diagonal of the matrix columnwise: ");
	// A = (double*)malloc(sizeof(double) * n * n);
	// for (int j = 0; j < n; j++) //iterate over columns
	// {
	// 	for (int i = j; i < n; i++) //iterate over rows
	// 	{
	// 		scanf("%lf", &A[j * n + i]);
  //     if(i!=j) A[i*n+ j] = A[j * n + i];
	// 	}
	// }

  n = 3;
  A = (double*)malloc(sizeof(double) * n * n);
  // A[0] = -4.;
  // A[1] =0.;
  // A[2] = 0.;
  // A[4] = 5.;
  // A[5] = 0.;
  // A[8] = -6.;
	A[0] = 1.;
	A[1] =2.;
	A[2] = -1.;
	A[4] = 1.;
	A[5] = 2.;
	A[8] = 1.;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      if(i != j) A[i*n+j] = A[j*n+i];
    }
  }

	printf("\nYour input matrix has the form: \n");
	printMatrix(n, n, A);
  printf("---------------------------------------------------------------\n");


	/*----------------------Tridiagonalization----------------------------------*/
  Q = (double*)malloc(sizeof(double) * n * n);
  /*set Q to the nxn unit matrix*/
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			Q[j * n + i] = 0;
			if (j == i) Q[j * n + i] = 1;
		}
	}

  double* w = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; i++) w[i] = Q[i];

  /*buffer Householder vector of current iteration step*/
  double* v = (double*)malloc(sizeof(double) * n);
  /*buffer column vector of A at current iteration step*/
  double* x = (double*)malloc(sizeof(double) * n);
  /*buffer row vector of A at current iteration step*/
  double* y = (double*)malloc(sizeof(double) * n);

  /*set x to the first column, y to the first row of A to start the iteration*/
  for (int i = 0; i < n; i++) v[i] = 0;
  for (int i = 0; i < n; i++) x[i] = A[i];

  /*tridiagonalize matrix A*/
  for (int j = 0; j < n-2; j++) //n-2 iteration steps
  {
    double abs = abs_vec(A, j+1, n, 0); //absolute value of part of column vector j of A

    /*calculate Householder vector v for iteration step j*/
		v[j+1] = A[j * n + (j+1)] + sgn(A[j * n + (j+1)]) * abs;
    for (int i = j+2; i < n; i++) v[i] = A[j * n + i];
		printVector(n,v);
		printf("\n");

    /*calculate P_j * A_j = A'_j*/
		double abs_v = abs_vec1(v, j+1, n, 1); //absolute value of HH vector v squared
    for (int l = j; l < n; l++) //iterate over columns
    {
      for (int k = j+1; k < n; k++) //iterate over rows
      {
        for (int d = j+1; d < n; d++)
        {
          A[l * n + k] += -2 / abs_v * v[k] * v[d] * x[d];
        }
      }
      /*set vector x to the next column of A*/
      for (int k = 0; k < n; k++) x[k] = A[(l + 1) * n + k];
    }

    printMatrix(n, n, A);

    /*calculate A'_j * P_j = P_j * A_j * P_j*/

    /*store first row vector of A' in y*/
    for (int l = 0; l < n; l++) y[l] = A[l*n];

    for (int l = j; l < n; l++) //iterate over columns
    {
      for (int k = j+1; k < n; k++) //iterate over rows
      {
        for (int d = j+1; d < n; d++)
        {
          A[k * n + l] += -2 /abs_v * v[k] * v[d] * y[d];
        }
      }
      /*set vector y to the next row of A*/
      for (int k = 0; k < n; k++) y[k] = A[k * n + (l + 1)];
    }

    printMatrix(n, n, A);
    printf("---------------------------------------------------------------\n");

    /*calculate matrix Q=P for each iteration step j*/
    for (int l = 0; l < n; l++) //iterate over columns of Q
    {
      //set vector w to transposed column vector
      for (int i = j; i < n; i++) w[i] = Q[i * n + l];
      /*calculate next Q-matrix analogeously to part 1 above, taking the transposition
      for the final result already into account*/
      for (int k = j; k < n; k++)
      {
        for (int d = j+1; d < n; d++) Q[k * n + l] += -2 / abs_v * v[k] * v[d] * w[d];
      }
    }

    /*reset vector x to the first column of the new matrix A*/
    for (int k = 0; k < n; k++) x[k] = A[(j + 2) * n + k];

    /*reset vector y to the first row of the new matrix A*/
    for (int l = 0; l < n; l++) y[l] = A[l * n + (j + 2)];

    /*reset HH-vector*/
    for (int i = 0; i < j+2; i++) v[i] = 0;

  }

  printf("Resulting tridiagonal matrix:\n");
  printMatrix(n, n, A);

  /*print orthogonal matrix Q (already transposed)*/
  printf("\nResulting orthogonal matrix Q: \n");
  printMatrix(n, n, Q);

	/*check result*/
	double* QA = (double*)malloc(sizeof(double) * n * n);
	double* QAQ = (double*)malloc(sizeof(double) * n * n);
	matrix_mul(Q, A, QA, n);
	matrix_mul(QA, Q, QAQ, n);
	printf("\nOriginal Matrix A: \n");
	printMatrix(n, n, QAQ);
	free(QA);
	free(QAQ);

	/*---------------------Calculate QR Decomposition of A----------------------*/
	// /*allocate nxn-matrix Q and set it to the unit matrix*/
	// Q = (double*)malloc(sizeof(double) * n * n);
	// for (int j = 0; j < n; j++)
	// {
	// 	for (int i = 0; i < n; i++)
	// 	{
	// 		Q[j * n + i] = 0;
	// 		if (j == i) Q[j * n + i] = 1;
	// 	}
	// }
  //
	// /*perform QR-decomposition*/
  // get_QR_decomposition(A, Q, n, m);
	/*--------------------------------------------------------------------------*/

	// /*print upper triangular matrix R with HH-vectors below its diagonal*/
	// printf("\nResulting matrix upper triangular matrix R with HH-vectors:\n");
	// printMatrix(n, m, A);
  //
	// /*print orthogonal matrix Q (already transposed)*/
	// printf("\nResulting orthogonal matrix Q: \n");
	// printMatrix(n, n, Q);

	/*free allocated memory*/

	free(A);
	free(Q);
  free(v);
  free(x);
  free(y);
  free(w);

	return 0;
}
