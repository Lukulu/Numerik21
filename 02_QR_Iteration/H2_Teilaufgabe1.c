/*@author: Laura Zywietz Rolon, Alina Becker*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <stdbool.h>

/*compile command:
gcc -Wall -pedantic -Wextra -Wno-unused-parameter H2_Teilaufgabe1.c -o ex -lm */

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

/*perform matrix multiplication of nxn matrices A, B and store
the result in AB*/
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

/*determine the sign of a given number a*/
int sgn(double a) {
	if (a >= 0) return 1;
	else return -1;
}

/*calculate minimum of two numbers a, b*/
int findMin(int a, int b) {
	if (a < b) return a;
	else return b;
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

    // /*write the normalized HH vector into thej-th column of A below the diagonal.*/
    // if (j != findMin(m - 1, n))
    // {
    //   for (int k = (j + 1); k < n; k++) A[j * n + k] = v[k] / v[j];
    // }
  }
  free(w);
  free(x);
  free(v);
}

/*calculate QR-decomposition of a band matrix A of dimension nxn with*/
void get_QR_decomp_band_matrix(double* A, double* Q, double* R, int n){
  /*allocate vector w of size n, and set it to the first column of Q*/
  double* w = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; i++) w[i] = Q[i];

  /*buffer column vector of A at current iteration step*/
  double* x = (double*)malloc(sizeof(double) * n);
  /*set x to the first column of A to start the iteration*/
  for (int i = 0; i < n; i++) x[i] = A[i];

  /*buffer Householder vector of current iteration step*/
  double* v = (double*)malloc(sizeof(double) * n);

	/*set R = A and Q = I to start the iteration*/
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			R[j*n+i] = A[j*n+i];
			Q[j * n + i] = 0;
			if (j == i) Q[j * n + i] = 1;
		}
	}

  /*perform QR-decomposition: calculate HH-vectors and R*/
  for (int j = 0; j < n; j++) //n iteration steps
  {
    double abs = abs_vec(R, j, n, 0); //absolute value of column vector j of A

    /*calculate Householder vector v for iteration step j*/
    for (int i = j; i < n; i++)
    {
      if (i == j) v[i] = R[j * n + i] + sgn(R[j * n + i]) * abs;
      else     v[i] = R[j * n + i];
    }

    //calculate H_j * A_j
    double abs_v = abs_vec1(v, j, n, 1); //absolute value of HH vector v squared
    for (int l = j; l < n; l++) //iterate over columns
    {
      for (int k = j; k < n; k++) //iterate over rows
      {
        for (int d = j; d < n; d++)
        {
          R[l * n + k] += -2 / abs_v * v[k] * v[d] * x[d];
        }
      }
      /*set vector x to the next column of A*/
      for (int k = 0; k < n; k++) x[k] = R[(l + 1) * n + k];
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
    for (int k = 0; k < n; k++) x[k] = R[(j + 1) * n + k];
  }
  free(w);
  free(x);
  free(v);
}

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
	int n = 0; //number of rows & columns

	/*----------------------Read Matrix from Command Line-----------------------*/
	/*Read number of rows n, number of columns m and the matrix A itself
	from console and store them in the corresponding containers.*/
	printf("Please enter the matrix size:\n");
	printf("Number of rows and columns: ");
	scanf("%d", &n);

	//terminate program if n = 0
	if (n == 0)
	{
		printf("Error! The number of rows and columns has to be different from zero.\n");
		return 1;
	}

	/*the matrix is columnwise read in to allow for a
	faster memory access of its columns.*/
	printf("Please enter the lower diagonal of the matrix columnwise: ");
	A = (double*)malloc(sizeof(double) * n * n);
	for (int j = 0; j < n; j++) //iterate over columns
	{
		for (int i = j; i < n; i++) //iterate over rows
		{
			scanf("%lf", &A[j * n + i]);
      if(i!=j) A[i*n+ j] = A[j * n + i];
		}
	}

	printf("\nYour input matrix has the form: \n");
	printMatrix(n, n, A);

	/*----------------------Tridiagonalization----------------------------------*/
	printf("-----------------------Tridiagonalization------------------------\n");
	double* P; //matrix P
	double eps = 1e-8; //check whether the HH-vector is zero, i.e. smaller than eps

  P = (double*)malloc(sizeof(double) * n * n); //store orthogonal transformation
  /*set P to the nxn unit matrix*/
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			P[j * n + i] = 0;
			if (j == i) P[j * n + i] = 1;
		}
	}

  double* w = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; i++) w[i] = P[i];

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
  for (int j = 0; j < n-2; j++) //n-1 iteration steps
  {
    double abs = abs_vec2(A, j, n, 0); //absolute value of part of column vector j of A

    /*calculate Householder vector v for iteration step j*/
		bool not_zero_entries = false;

    for (int i = j+1; i < n; i++)
		{
			 if(i == (j+1)) v[i] = A[j * n + i] + sgn(A[j * n + i]) * abs;
			 else v[i] = A[j * n + i];
			 if(fabs(v[i]) > eps) not_zero_entries = true;
		}

		/*terminate iteration in the case that the HH vector is zero
		-> Input matrix is already diagonal*/
		if(!not_zero_entries)
		{
			printf("Matrix already diagonal. HH-vector becomes zero. Iteration is terminated.\n");
			break;
		}

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

    /*calculate A'_j * P_j = P_j * A_j * P_j*/
    /*store first row vector of A' in y*/
    if(j==0) for (int l = 0; l < n; l++) y[l] = A[l*n];

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

    /*calculate matrix P for each iteration step j*/
    for (int l = j; l < n; l++) //iterate over columns of P
    {
      for (int k = j+1; k < n; k++)
      {
        for (int d = j+1; d < n; d++) P[l * n + k] += -2 / abs_v * v[k] * v[d] * w[d];
      }
			for (int k = 0; k < n; k++) w[k] = P[(l+1) * n + k];
    }

    /*reset vector x to the first column of the new matrix A*/
    for (int k = 0; k < n; k++) x[k] = A[(j + 1) * n + k];

    /*reset vector y to the first row of the new matrix A*/
    for (int l = 0; l < n; l++) y[l] = A[l * n + (j + 1)];

		/*reset vector w to the first column of P*/
		for (int k = 0; k < n; k++) w[k] = P[(j + 1) * n + k];

    /*reset HH-vector*/
    for (int i = 0; i < j+2; i++) v[i] = 0;

  }

	/*print tridiagonalized matrix A*/
  printf("\nResulting tridiagonal matrix A1:\n");
  printMatrix(n, n, A);

  /*print orthogonal matrix P*/
  printf("Resulting orthogonal matrix P:\n");
  printMatrix(n, n, P);

	/*check result by calculating P*A*P to get the original matrix back*/
	double* Pt = (double*)malloc(sizeof(double) * n * n);
	transposeMatrix(n, P, Pt);
	double* PtA = (double*)malloc(sizeof(double) * n * n);
	double* PtAP = (double*)malloc(sizeof(double) * n * n);
	matrix_mul(Pt, A, PtA, n);
	matrix_mul(PtA, P, PtAP, n);
	printf("Original Matrix A has the form: \n");
	printMatrix(n, n, PtAP);
	free(PtA);
	free(Pt);
	free(PtAP);
	/*--------------------------------------------------------------------------*/

	/*------------------------------QR-Iteration--------------------------------*/
	bool no_zero_entry = false;
	double epsilon = 0.00000001;
	/*allocate nxn-matrix Qk for iteration step k*/
	double* Qk = (double*)malloc(sizeof(double) * n * n);
	/*allocate nxn-matrix Rk for iteration step k*/
	double* Rk = (double*)malloc(sizeof(double) * n * n);
	/*allocate nxn-matrix Q0 for iteration step 0 and set it to the unit matrix*/
	double* Q0 = (double*)malloc(sizeof(double) * n * n);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			Q0[j * n + i] = 0;
			if (j == i) Q0[j * n + i] = 1;
		}
	}
	/*allocate nxn-matrix Q and set it to the zero matrix*/
	double* Q = (double*)malloc(sizeof(double) * n * n);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++) Q[j * n + i] = 0;
	}

	int k = 0;
	while (k >= 0)
	{
		get_QR_decomp_band_matrix(A, Qk, Rk, n);
		matrix_mul(Rk, Qk, A, n);

		/*calculate matrix Q=Q_1*..*Q_k for iteration step k*/
		// matrix_mul(Qk, Q0, Q, n);
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++) Q[j * n + i] = 0;
		}

		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
				for (int d = 0; d < n; d++) Q[j*n+i] += Qk[d*n+i]*Q0[j*n+d];
				/*check whether entries of Qk next to the diagonal
				are smaller than eps=1e-8*/
				if(i != j)
				{
					if (fabs(Qk[j*n+i]) > epsilon) no_zero_entry = true;
				}
			}
		}
		/*terminate iteration if all entries of Qk apart from the diagonal are
		smaller than eps = 1e-8*/
		if (!no_zero_entry) break;

		/*copy Q into Q0 to start the next iteration step*/
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				Q0[j*n+i] = Q[j*n+i];
			}
		}
		k++; //increase iterations step
		no_zero_entry = false;
	}

	/*print result QR iteration*/
	printf("--------------------------QR-Iteration---------------------------\n");
	printf("\nResulting matrix QR-Iteration R:\n");
	printMatrix(n, n, A);

	// double* Qt = (double*)malloc(sizeof(double) * n * n);
	// transposeMatrix(n, Q, Qt);
	printf("\nResulting matrix Q:\n"); //->transposed??
	printMatrix(n, n, Q);
	/*--------------------------------------------------------------------------*/

	/*--------------------------Calculate Eigenvectors--------------------------*/
	printf("-----------------Eigenvectors & Eigenvalues----------------------\n");
	/*calculate eigenvectors from P*Q from part 1,2 and store the result in Q0*/
	matrix_mul(P, Q, Q0, n);

	/*print P*Q*/
	printf("\nResulting matrix PQ:\n");
	printMatrix(n, n, Q0);

	/*print eigenvectors (from columns of PQ) and corresponding eigenvalues
	(from the diagonal of A) obtained from QR-iteration*/
	for (int j = 0; j < n; j++)
	{
		printf("%d. Eigenvector to eigenvalue %4.2lf :", j, A[j*n+j]);
		for (int i = 0; i < n; i++)
		{
			if(i == 0) printf("( ");
			if(i==(n-1)) printf(" %5.3lf )\n", Q0[j*n+i]);
			else printf(" %5.3lf ,", Q0[j*n+i]);
		}
		printf("\n");
	}
	printf("-----------------------------------------------------------------\n");
	/*--------------------------------------------------------------------------*/

	/*free allocated memory*/
	free(Qk);
	free(Rk);
	free(Q0);
	free(Q);
	free(A);
	free(P);
  free(v);
  free(x);
  free(y);
  free(w);

	return 0;
}
