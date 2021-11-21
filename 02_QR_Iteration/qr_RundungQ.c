# include <stdio.h>
# include <stdlib.h>
# include <math.h>


/*-----------------------Funktionen-----------------------*/
/*Matrix A ausgeben*/
void printMatrix(int n, double* A) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j == 0) printf("( ");
			printf("%8.8lf ", A[j * n + i]);
			if (j == (n - 1)) printf(")\n");
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

/*determine the sign of a given number a*/
int sgn(double a) {
	if (a >= 0) return 1;
	else return -1;
}

double** calculate_Q(double* A, int n, double* Q, double* v, double* w, double* x, double** QR){
  //QR[0] = Q, QR[1] = R

	//Q auf Einheitsmatrix setzen
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			Q[j * n + i] = 0;
			if (j == i) Q[j * n + i] = 1;
		}
	}

	for (int i = 0; i < n; i++) w[i] = Q[i];

	//set x to the first column of A to start the iteration
	for (int i = 0; i < n; i++) x[i] = A[i];

	/*perform QR-decomposition: calculate HH-vectors and R*/
	for (int j = 0; j <= n; j++) //min(m-1, n)+1 iteration steps
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
		for (int l = j; l < n; l++) //iterate over columns
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

		/*reset vector x to the first column of
		the new matrix A*/
		for (int k = 0; k < n; k++) x[k] = A[(j + 1) * n + k];

		/*write the normalized HH vector into the
		j-th column of A below the diagonal.*/
		if (j != n)
		{
			for (int k = (j + 1); k < n; k++) A[j * n + k] = v[k] / v[j];
		}

		/*for debugging*/
		// printf("Iterationstep: %d\n", j);
		//
		// printf("\nResulting matrix of the iteration step %d with HH vectors: \n", j);
		// printMatrix(n, m, A);
		//
		// printf("\nResulting matrix Q of the iteration step %d: \n", j);
		// printMatrix(n, n, Q);
	}

	/*--------------------------------------------------------------------------*/
	/*print upper triangular matrix R with HH-vectors below its diagonal*/
	printf("\nResulting matrix upper triangular matrix R with HH-vectors:\n");
	printMatrix(n, A);
	QR[1] = A;

	/*print orthogonal matrix Q (already transposed)*/
	printf("\nResulting orthogonal matrix Q: \n");
	printMatrix(n, Q);
	QR[0] = Q;

  return QR;
}

double* calculate_R(double* A, int n, double* v, double* x){

	//set x to the first column of A to start the iteration
	for (int i = 0; i < n; i++) x[i] = A[i];

	/*perform QR-decomposition: calculate HH-vectors and R*/
	for (int j = 0; j <= n; j++) //min(m-1, n)+1 iteration steps
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
		for (int l = j; l < n; l++) //iterate over columns
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

		/*reset vector x to the first column of
		the new matrix A*/
		for (int k = 0; k < n; k++) x[k] = A[(j + 1) * n + k];

		/*write the normalized HH vector into the
		j-th column of A below the diagonal.*/
		if (j != n)
		{
			for (int k = (j + 1); k < n; k++) A[j * n + k] = v[k] / v[j];
		}

		/*for debugging*/
		// printf("Iterationstep: %d\n", j);
		//
		// printf("\nResulting matrix of the iteration step %d with HH vectors: \n", j);
		// printMatrix(n, m, A);
	}

	/*--------------------------------------------------------------------------*/
	/*print upper triangular matrix R with HH-vectors below its diagonal*/
	printf("\nResulting matrix upper triangular matrix R with HH-vectors:\n");
	printMatrix(n, A);

  return A;
}

/*-----------------------Main-Methode-----------------------*/
int main(int argc, char** argv){
  double* A; //Eingabematrix A (symmetrisch)
  int n = 0;
	double* Q;
	double* R;
	double** QR;

	/*für QR Zerlegung benötigt:*/
	double* v; //HH vector of current iteration step j
	double* x; //j-th column of A
	double* w; //store 1,...,m unit vectors of R^n temporarily

  QR = (double**)malloc(sizeof(double*)*2);
  A = (double*)malloc(sizeof(double)*n*n);
	Q = (double*)malloc(sizeof(double)*n*n);
	R = (double*)malloc(sizeof(double)*n*n);
	v = (double*)malloc(sizeof(double)*n); //buffer Householder vector of current iteration step
	x = (double*)malloc(sizeof(double)*n); //buffer column vector of A at current iteration step
	w = (double*)malloc(sizeof(double)*n);

  printf("Bitte die Matrixdimension eingeben:\n");
	scanf("%d", &n);

  //terminate program if n = 0
	if (n == 0)
	{
		printf("Fehler: die Dimension darf nicht null sein.\n");
		return 1;
	}

  //Matrix spaltenweise einlesen (linke untere Dreiecksmatrix)
  printf("Matrix ab 1.Zeile bis Hauptdiagonale spaltenweise eingeben:\n");
  for (int i = 0; i < n; i++){
    for (int j = 0; j <= i; j++){
      scanf("%lf", &A[i * n + j]); //Spalten über (einschließlich) Hauptdiagonalen
      if (i != j){
        A[j * n + i] = A[i * n + j]; //Spalten unter Hauptdiagonalen
      }
    }
  }

  printf("\nEingabematrix sieht folgendermaßen aus:\n");
	printMatrix(n, A);

  for (int c = 0; c < 1; c++){
	//for (int c = 0; c < n; c++){ //öfters rechnen, damit A sich immer mehr triangularisiert, am Ende stehen EW auf der Diagonalen
    printf("Q und R werden berechnet:\n");
		QR = calculate_Q(A, n, Q, v, w, x, QR);
		R = QR[1];
		Q = QR[0];
		printf("Q und R wurden berechnet.\n");

		/*linkes unteres Dreieck von R auf 0 setzen*/
		for (int i = 0; i < n; i++){
	    for (int j = 0; j < i; j++){
	      if (i != j){
	        R[j * n + i] = 0; //Spalten unter Hauptdiagonalen
	      }
	    }
	  }
		printf("\nMatrix R ohne HH Vektoren:\n");
		printMatrix(n, R);

		printf("Matrix Q:\n");
		printMatrix(n, Q);

    printf("Matrix A(k) auf Null setzen, um A(k+1) zu berechnen:\n");
		for (int i = 0; i < n; i++){
	    for (int j = 0; j <= i; j++){
	      A[i * n + j] = 0; //Spalten über (einschließlich) Hauptdiagonalen
	      if (i != j){
	        A[j * n + i] = A[i * n + j]; //Spalten unter Hauptdiagonalen
	      }
	    }
	  }

    printf("A(k+1) = R(k)*Q(k) wird berechnet...\n");
    /*A(k+1) = R(k)*Q(k) berechnen: */
		for (int j = 0; j < n; j++){
			for (int i = 0; i < n; i++){
				for (int k = 0; k < n; k++){
					A[i * n + j] = A[i * n + j] + (R[i * n + k]*Q[k * n + j]);
				}
			}
		}
		printf("Matrix A(k+1):\n");
		printMatrix(n, A);
	}

  //allokierten Speicher freigeben
  free(A);
	printf("A wurde freigegeben\n");
	free(Q);
	printf("Q wurde freigegeben\n");
	free(R);
	printf("R wurde freigegeben\n");
	free(QR);
	printf("QR wurde freigegeben\n");
	free(v);
	printf("v wurde freigegeben\n");
	free(x);
	printf("x wurde freigegeben\n");
	free(w);
	printf("w wurde freigegeben\n");
	return 0;
}
