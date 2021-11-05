# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/*Print the matrix for debugging purpose*/
void printMatrix(int n, int m, double *A){
	printf("\nThe matrix has the form:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			printf("%lf ", A[j*n+i]);
			if (j==(m-1)){ printf("\n");}
		}
	}
	printf("\n");
}

void printVector(int length, double*v){
	printf("The vector has the form:\n");
	for (int i = 0; i < length; i++) {
		if(i==0) printf("( ");
		printf("%lf ", v[i]);
		if(i==(length-1)) printf(")\n");
	}
}

int findMin(int a, int b){
	if (a < b) return a;
	else return b;
}

int sgn(double a){
	if(a >= 0) return 1;
	else return -1;
}

double abs_vec(double *A, int pos, int length, int squared){
	double res = 0;
	for (int i = pos; i < length; i++) {
		res += A[pos*length+i]*A[pos*length+i];
	}
	if(squared) return res;
	else return sqrt(res);
}

double abs_vec1(double *v, int pos, int length, int squared){
	double res = 0;
	for (int i = pos; i < length; i++) {
		res += v[i]*v[i];
	}
	if(squared) return res;
	else return sqrt(res);
}


int main (int argc, char** argv)
{
	/*Read number of rows n, number of columns m and the matrix A itself
	from console and store them in the corresponding containers.*/
	double *A;
	double *v;
	double *x;
	int n = 0, m = 0;

	// printf("Please enter the matrix size:\n");
	// printf("Number of rows: ");
	// scanf("%d", &n);
	// printf("Number of columns: ");
	// scanf("%d", &m);
	// if(n == 0 || m == 0){
	// 	printf("Error! The number of rows and columns has to be different from zero. Try again.\n");
	// 	return 1;
	// }
	//
	// printf("Please enter the matrix column-wise:\n ");
	// A = (double *) malloc(sizeof(double) * n* m);
	// for (int j = 0; j < m; j++) {
	// 	for (int i = 0; i < n; i++) {
	// 		scanf("%lf", &A[j*n+i]);
	// 	}
	// }
	n = 3;
	m = 3;
	A = (double *) malloc(sizeof(double) * n* m);
	A[0] = 2;
	A[1] = 1;
	A[2] = 2;
	A[3] = 4;
	A[4] = 1;
	A[5] = -3;
	A[6] = -4;
	A[7] = 2;
	A[8] = 0;

	printMatrix(n, m, A);


	/*Perform QR-decomposition*/
	v = (double *) malloc(sizeof(double)*n);
	x = (double *) malloc(sizeof(double)*n);
	for (int i = 0; i < n; i++) {
		x[i] = A[i];
	}

	for (int j = 0; j < findMin(m-1, n); j++)
	{
		printf("Iterationstep: %d\n", j);
		printVector(n, x);

		double abs = abs_vec(A,j,n, 0);
		printf("Betrag von a: %2.1lf \n", abs);
		//calculate Householder vector v for iteration step j
		for (int i = j; i < n; i++)
		{
			printf("A[1,i]: %lf\n", A[j*n+i]);
			if(i==j) v[i] = A[j*n+i] + sgn(A[j*n+i])*abs;
			else     v[i] = A[j*n+i];//v[j];
		}
		//v[j] = 1.0;
		printVector(n, v);

		//calculate R
		double abs_v = abs_vec1(v,j,n, 1);
		printf("Betrag v: %lf\n", abs_v);
		for (int l = j; l < m; l++) { //columns
			for (int k = j; k < n; k++) { //rows
				for (int d = j; d < n; d++) {
					// printf("%lf\n", A[l*n+k]);
					if(k==d) A[l*n+k] += -2/abs_v*v[k]*v[k]*x[k];
					else     A[l*n+k] += -2/abs_v*v[k]*v[d]*x[d];
				}
			}
			for (int k = 0; k < n; k++) {
				x[k] = A[(l+1)*n+k];
			}
		}
		for (int k = 0; k < n; k++) {
			x[k] = A[(j+1)*n+k];
		}

		if(j!=(m-1)){
			for (int k = (j+1); k < n; k++) {
				A[j*n+k] = v[k]/v[j];
			}
		}
		printMatrix(n, m, A);
	}
	/*free allocated memory*/
	free(A);
	free(v);
	free(x);
	return 0;
}
