# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/*Matrix A ausgeben*/
void printMatrix(int n, double* A) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j == 0) printf("( ");
			printf("%3.3lf ", A[j * n + i]);
			if (j == (n - 1)) printf(")\n");
		}
	}
	printf("\n");
}

int main(int argc, char** argv){
  double* A; //Eingabematrix A (symmetrisch)
  int n = 0;

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
  A = (double *) malloc(sizeof(double)*n*n);
  for (int i = 0; i < n; i++){
    for (int j = 0; j <= i; j++){
      scanf("%lf", &A[i * n + j]); //Spalten über Hauptdiagonalen
      if (i != j){
        A[j * n + i] = A[i * n + j]; //Spalten unter Hauptdiagonalen
      }
    }
  }

  printf("\nEingabematrix sieht folgendermaßen aus:\n");
	printMatrix(n, A);

/*
  double * R, Q;
  Q = calculate_Q(A, n);
  R = calculate_R(A, n);
*/

  //allokierten Speicher freigeben
  free(A);
}
