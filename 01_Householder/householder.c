//Header-Dateien
#include <stdio.h>
#include <stdlib.h>
#include <string.h> //für strtok
#include <math.h>
/*
int getSign(int sign)
{
    if (sign < 0)
    {
        return -1;
    }
    else if (sign >= 0)
    {
        return +1;
    }
}

double getNorm(double* x)
{
    int sqrt_norm = 0;
    double norm = 0;
    int x_lentgh = sizeof(x) / sizeof(double);
    for (int l = 0; l < x_lentgh; ++l)
    {
        sqrt_norm = sqrt_norm + pow(x[l], 2.);
        double norm = pow(sqrt_norm, 0.5);
    }
    return norm;
}
*/
int main (int argc, char **argv){
    // (mxn)-Matrix einlesen:
    // 1) m und n bestimmen
    double* dim;
    dim = (double *) malloc(2*sizeof(double));
    printf("m = ");
    scanf("%lf", &dim[0]);
    printf("n = ");
    scanf("%lf", &dim[1]);
    printf("m ist %lf und n ist %lf\n", dim[0], dim[1]);

    if (dim[0] < dim[1])
    {
        printf("n darf nicht größer sein als m!\n");
        free(dim);
    }

    // 2) Matrixeinträge einlesen (_ab)
    int wert_anzahl = dim[0] * dim[1];
    double* a;
    double** b;
    b = (double**) malloc(dim[0] * sizeof(double*));
    a = (double*) malloc(dim[1] * sizeof(double));
    char zahlen_eingabe[8];
    char zahlen_matrix[10];
    printf("%d Matrixeinträge mit Leerzeichen getrennt eingeben (wird zeilenweise in Matrix geschrieben): \n", wert_anzahl);
    /*scanf("%c\n", zahlen_eingabe);
    zahlen_matrix = strtok(zahlen_eingabe, " ");
    printf("%lf\n", (double) zahlen_matrix[0]);
*/
    double matrix[(int) dim[0]][(int) dim[1]];
    for (int i = 0; i < dim[0]; ++i)
    {
        for (int j = 0; j < dim[1]; ++j)
        {
            scanf("%lf", &matrix[i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i < dim[0]; ++i)
    {
        for (int j = 0; j < dim[1]; ++j)
        {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
/*
    int count = 0;
    //for (int i = 0; i < wert_anzahl; ++i)
    while (count <= wert_anzahl-1)
    {
        printf("%d.: ", count);
        scanf("%s", &zahlen_matrix[count]);
        printf("\n");
        count += 1;
    }
    printf("%s\n", zahlen_matrix);
    */

    /*
    for (int i = 1; i < wert_anzahl; ++i)
    {
        zahlen_matrix[i] = (double) strtok(NULL, " "); // delimiter=" "
    }
    printf("Zahlen_matrix %lf, %lf, %lf, %lf\n", zahlen_matrix[0], zahlen_matrix[1], zahlen_matrix[2], zahlen_matrix[3]); //Sollte erste eingegeben Zahl printen
    */
    /*while (zahlen_matrix != NULL)
    {
        printf("%s\n", zahlen_matrix);
        zahlen_matrix = strtok(NULL, " "); // TODO
    }*/
/*
    // in zahlen_matrix stehen nun die Matrixeinträge der zu zerlegenden Matrix A
    for (int i = 0; i < dim[0]; ++i)
    {
        for (int j = 0; j < dim[1]; ++j)
        {
            // a wird für jede Zeile der Matrix neu überschrieben
            int m = dim[0];
            a[j] = zahlen_matrix[i*m + j];
            printf("zahlen_matrix: %d\n", zahlen_matrix[i*m + j]);
            printf("a = %f\n", a[j]);
        }
        // a ist n lang (#Spalten) und wird in b gespeichert. Das passiert m mal (#Zeilen).
        b[i] = a;
    }
    for (int i = 0; i < dim[0]; ++i)
    {
        printf("b = %lf\n", *b[i]);
    }
    */

    // Berechne R

/*
    for (int k = 1; k <= dim[1]; ++k)
    {
        int counter = 0;
        double x[(int) dim[0]-k];
        for (int t = k; t <= dim[0]; ++t)
        {
            x[counter] = b[t* (int) dim[0]][k]; //ab der k-ten Zeile (l=k) bis zur m-ten Zeile immer den Wert der k-ten Spalte auswählen
            counter += 1;
        }
        int sign = getSign(x[0]);
        int norm = getNorm(x);
        //Bestimme Dimension Einheitsvektor
        double e_1[(int) dim[0]-k];
        for (int i = 0; i < k; ++i)
        {
            if (i == 0)
            {
                e_1[i] = 1;
            } else if (i != 0)
            {
                e_1[i] = 0;
            }
        }
        double term = e_1[0] * sign * norm;
        x[0] = x[0] + term;

        double* v = x;
    }
*/
    free(dim);
    free(b);
    free(a);
    return 0;
/*    for (int i = 1; i < argc; ++i) //erstes Argument ist ./a.out, aber das soll ignoriert werden
    {
        printf("argv[%d]: %s\n", i, argv[i]);
    }
*/
}
