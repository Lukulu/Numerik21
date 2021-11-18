# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/*---------------functions---------------*/

//print the matrix A
void printMatrix(int m, int n, double *A){
	
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf("%3.3lf ", A[j*m+i]);
			if (j==(n-1)){ printf("\n");}
		}
	}
	printf("\n");
}


//print vector v

void printVector(int length, double*v){
	printf("The vector is:\n");
	for (int i = 0; i < length; i++) {
		if(i==0) printf("( ");
		printf("%3.3lf ", v[i]);
		if(i==(length-1)) printf(")\n");
	}
}
/*minimum*/
int findMin(int a, int b){
	if (a < b) return a;
	else return b;
}

/*sign of a given number a*/
int sgn(double a){
	if(a >= 0) return 1;
	else return -1;
}

/*determine absolute value of column vector at pos of A
of length n according to the Euclidean norm.
If squared = 1: return the squared value to avoid
its calculation afterwars by squaring the root.*/
double abs_vec(double *A, int pos, int m, int squared){
	double res = 0;
	for (int i = pos; i < m; i++) {
		res += A[pos*m+i]*A[pos*m+i];
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

/*SKP*/
double SKP(double*v, double*Q, int m, int n, int j, int pos){
	double res = 0;
	for (int i = j; i < m; i++) {
		res += v[i]*Q[pos*m+i];
	}
	return res;                                                        
}
/*---------------------------------------*/

int main (int argc, char** argv)
{
	double *A; //input A
	double *v; //HH vector
	double *x; 
	double *y;
	double *Q;
	int n = 0, m = 0;

	//Input
	printf("Please enter the matrix size:\n");
	printf("Number of rows: ");
	scanf("%d", &m);
	printf("Number of columns: ");
	scanf("%d", &n);
	if(n == 0 || m == 0)
	{
		printf("Error! The number of rows and columns has to be different from zero. Try again.\n");
		return 1;
	}

	printf("Please enter the matrix column-wise:");
	A = (double *) malloc(sizeof(double) * n* m);
	for (int j = 0; j < n; j++) 
	{
		for (int i = 0; i < m; i++) 
		{
			scanf("%lf", &A[j*m+i]);
		}
	}

	printf("Your input matrix has the form: \n");
	printMatrix(m, n, A);

//Q as 1
	Q = (double *) malloc(sizeof(double)*m*m);
	for (int j = 0; j < m; j++) {
		for (int i = 0; i < m; i++) {
			if(i==j) Q[j*m+i] = 1;
			else Q[j*m+i] = 0;
		}
	}

	v = (double *) malloc(sizeof(double)*m);
	x = (double *) malloc(sizeof(double)*m);
	y = (double *) malloc(sizeof(double)*m);
    for (int i = 0; i < n; i++) x[i] = A[i];
    for (int i = 0; i < n; i++) y[i] = Q[i*m];
    
	/*perform QR-decomposition:
	calculate HH-vectors and R*/
	for (int j = 0; j <= findMin(n-1, m); j++)
	{
		printf("Iterationstep: %d\n", j);

		double abs = abs_vec(A,j,m, 0);                                          //norm of jth column
		printf("Absoulte value of column vector j of A: %2.1lf \n", abs);

		/* Householder vector v for iteration j*/
		for (int i = j; i < m; i++)
		{
			if(i==j) v[i] = A[j*m+i] + sgn(A[j*m+i])*abs;
			else     v[i] = A[j*m+i];
		}
		printVector(n, v);

		
		double abs_v = abs_vec1(v,j,m, 1); //absolute value of HH vector v squared
		printf("Absolute value HH vector v: %2.1lf\n", abs_v);
		//calculate H_j * A_j
		for (int l = j; l < n; l++) //iterate over columns
		{			
			for (int k = j; k < m; k++)
			{
			    for (int d = j; d < m; d++) 
			    {
			        if(k==d) A[l*m+k] += -2/abs_v*v[k]*v[k]*x[k];
				    else     A[l*m+k] += -2/abs_v*v[k]*v[d]*x[d];
				    
				}
			}
			for (int k = 0; k < m; k++) x[k] = A[(l+1)*m+k];
		}
        //calculate Q
        for(int l = 0; l<n ;l++){
			for (int k = j; k < m; k++) y[k] = Q[k*m+l];
            for (int d = j; d < m; d++)
			{
			    for (int k = j; k < m; k++)
			    {
			        Q[m*d+l] += -2/abs_v*v[k]*v[d]*y[k];
				}
			}
        }
		for (int k = 0; k < n; k++) x[k] = A[(j+1)*n+k];
		
		/*writing the normalized HH vector in A.*/
		if(j!=(m-1))
		{
			for (int k = (j+1); k < m; k++) A[j*m+k] = v[k]/v[j];
		}

		printf("\nResulting matrix of the iteration step %d with HH vectors: \n", j);
		printMatrix(m, n, A);

		printf("\nQ in  %d iteration \n",j );
		printMatrix(m, n, Q);
	}
	free(A);
	free(v);
	free(x);
	free(Q);
    free(y);
	return 0;
	}