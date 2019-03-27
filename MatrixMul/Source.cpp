#define _CRT_SECURE_NO_WARNINGS

// Celem tego programu jest prezentacja pomiaru i analizy 
//efektywnosci programu za pomoc¹  CodeAnalyst(tm).
// Implementacja mno¿enia macierzy jest realizowana za pomoca typowego 
// algorytmu podrêcznikowego. 
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include "omp.h"

#define USE_MULTIPLE_THREADS true
#define MAXTHREADS 128
int NumThreads;
double start;

static const int ROWS = 1000;     // liczba wierszy macierzy
static const int COLUMNS = 1000;  // lizba kolumn macierzy

float matrix_a[ROWS][COLUMNS];    // lewy operand 
float matrix_b[ROWS][COLUMNS];    // prawy operand
float matrix_r[ROWS][COLUMNS];    // wynik

FILE *result_file;

void initialize_matrices()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
//#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_a[i][j] = (float)rand() / RAND_MAX;
			matrix_b[i][j] = (float)rand() / RAND_MAX;
			matrix_r[i][j] = 0.0;
		}
	}
}

void initialize_matricesZ()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_r[i][j] = 0.0;
		}
	}
}
void print_result()
{
	// wydruk wyniku
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			fprintf(result_file, "%6.4f ", matrix_r[i][j]);
		}
		fprintf(result_file, "\n");
	}
}

void multiply_matrices_IKJ_sequence()
{
	// mnozenie macierzy 
//#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++)
		for (int k = 0; k < COLUMNS; k++)
			for (int j = 0; j < COLUMNS; j++)
				matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];

}

void multiply_matrices_JIK_sequence()
{
	for (int j = 0; j < COLUMNS; j++) {
		for (int i = 0; i < ROWS; i++) {
			float sum = 0.0;
			for (int k = 0; k < COLUMNS; k++) {
				sum = sum + matrix_a[i][k] * matrix_b[k][j];
			}
			matrix_r[i][j] = sum;
		}
	}
}

void multiply_matrices_JIK()
{
#pragma omp parallel for
	for (int j = 0; j < COLUMNS; j++) {
		for (int i = 0; i < ROWS; i++) {
			float sum = 0.0;
			for (int k = 0; k < COLUMNS; k++) {
				sum = sum + matrix_a[i][k] * matrix_b[k][j];
			}
			matrix_r[i][j] = sum;
		}
	}
}

void multiply_matrices_six_loops_IJK_II_KK_JJ_sequence() {
	int r = 25;
	for (int i = 0; i < ROWS; i += r)
		for (int j = 0; j < COLUMNS; j += r)
			for (int k = 0; k < COLUMNS; k += r)// kolejne fragmenty  
				for (int ii = i; ii < i+r; ii++)//fragment wyniku 
					for (int kk = k; kk < k+r; kk++) 
						for (int jj = j; jj < j+r; jj++) 
							matrix_r[ii][jj] += matrix_a[ii][kk] * matrix_b[kk][jj];
}

void multiply_matrices_six_loops_IJK_II_KK_JJ() {
	int r = 25;
#pragma omp parallel for
	for (int i = 0; i < ROWS; i += r)
		for (int j = 0; j < COLUMNS; j += r)
			for (int k = 0; k < COLUMNS; k += r)// kolejne fragmenty  
				for (int ii = i; ii < i + r; ii++)//fragment wyniku 
					for (int kk = k; kk < k + r; kk++)
						for (int jj = j; jj < j + r; jj++)
							matrix_r[ii][jj] += matrix_a[ii][kk] * matrix_b[kk][jj];
}

void print_elapsed_time()
{
	double elapsed;
	double resolution;

	// wyznaczenie i zapisanie czasu przetwarzania
	elapsed = (double)clock() / CLK_TCK;
	resolution = 1.0 / CLK_TCK;
	printf("Czas: %8.4f sec \n",
		elapsed - start);

	fprintf(result_file,
		"Czas wykonania programu: %8.4f sec (%6.4f sec rozdzielczosc pomiaru)\n",
		elapsed - start, resolution);
}

int main(int argc, char* argv[])
{
	//	 start = (double) clock() / CLK_TCK ;
	if ((result_file = fopen("classic.txt", "a")) == NULL) {
		fprintf(stderr, "nie mozna otworzyc pliku wyniku \n");
		perror("classic");
		return(EXIT_FAILURE);
	}


	//Determine the number of threads to use
	if (USE_MULTIPLE_THREADS) {
		SYSTEM_INFO SysInfo;
		GetSystemInfo(&SysInfo);
		NumThreads = SysInfo.dwNumberOfProcessors;
		if (NumThreads > MAXTHREADS)
			NumThreads = MAXTHREADS;
	}
	else
		NumThreads = 1;
	fprintf(result_file, "Klasyczny algorytm mnozenia macierzy, liczba watkow %d \n", NumThreads);
	printf("liczba watkow  = %d\n\n", NumThreads);

	initialize_matrices();
	start = (double)clock() / CLK_TCK;
	//multiply_matrices_JIK();
	//printf("JIK ");
	multiply_matrices_six_loops_IJK_II_KK_JJ_sequence();
	printf("IJK II KK JJ ");
	print_elapsed_time();

	fclose(result_file);

	return(0);
}
