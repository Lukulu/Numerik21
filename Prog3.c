# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <stdbool.h>

/*--------------------------------functions-----------------------------------*/
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

int sgn(double a) {
	if (a >= 0) return 1;
	else return -1;
}

int findMin(int a, int b) {
	if (a < b) return a;
	else return b;
}

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

double abs_vec1(double* v, int pos, int length, int squared) {
	double res = 0;
	for (int i = pos; i < length; i++) {
		res += v[i] * v[i];
	}
	if (squared) return res;
	else return sqrt(res);
}
/*----------------------------------------------------------------------------*/



int main(int argc, char** argv)
{
	double* A;                                        
	int n = 0;                                        
	printf("Gebe eine s.p.d. Matrix ein:\n");
	printf("Matrixgröße: ");
	scanf("%d", &n);
	printf("Einträge spaltenweise eingeben: ");
	A = (double*)malloc(sizeof(double) * n * n);
	for (int j = 0; j < n; j++)
	{
		for (int i = j; i < n; i++)
		{
			scanf("%lf", &A[j * n + i]);
      if(i!=j) A[i*n+ j] = A[j * n + i];
		}
	}
	printf("\nYour input matrix has the form: \n");
	printMatrix(n, n, A);
	double* b;                                         //input b
	double* x;                                         //input x
	int Eta;
	printf("b eingeben: ");
	b = (double*)malloc(sizeof(double) *n);
	for (int j = 0; j < n; j++){
	    scanf("%lf", &b[j]);
	}
	printf("x eingeben: ");
	x = (double*)malloc(sizeof(double) *n);
	for (int j = 0; j < n; j++){
	    scanf("%lf", &x[j]);
	}
	printMatrix(1, n, b);
	printMatrix(1, n, x);
//Alle Inputs


	double* r;
	r = (double*)malloc(sizeof(double) * n);
	double* Ax;
	Ax= (double*)malloc(sizeof(double) * n);
	matrix_mul(A, x, Ax, n);
	for (int j=0; j<n;j++){
	    r[j]=b[j]-Ax[j];
	}
	printMatrix(1, n, r);                                  //r0
	
	//while Norm r< eta*Norm b : rk+1=rk-akApk(aus Tutorium)

	/*free allocated memory*/
	free(A);
	free(b);
	free(x);
	free(r);
	free(Ax);
}