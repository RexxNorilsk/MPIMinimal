#include <iostream>
#include "mpi.h"
#include <Windows.h>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <fstream>
//Вариант 9
// 
//mpiexec -n 4 "C:\Parallel 3.2\MinimalGrid\x64\Debug\MinimalGrid.exe"
//mpiexec -n 4 "B:\3.1 Parralel\GridMinimal\MPIMinimal\x64\Debug\MinimalGrid.exe"
// 
//Топология
#define ndim 1

//Шаги
#define hEnd 0.000025
#define hStart 0.01

//Область
#define top 5
#define right 5
#define left -5
#define bottom -5


using namespace std;

struct point {
	double x;
	double y;
	double value;
};

int randomRange(int max, int min = 1) {
	return min + rand() % (max - min);
}

double func(point data) {
	return 9.2* data.x - 0.6 * data.y + exp(0.81 * data.x * data.x + 0.19 * data.y * data.y);
	//9.2*x-0.6*y+exp(0.81*x^2+0.19*y^2)
}

point findMin(double startElementX, double endElementX, double startElementY, double endElementY, double step) {
	point min, temp;
	min.x = startElementX * step;
	min.y = startElementY * step;
	min.value = func(min);
	for (double i = startElementX; i < startElementX + endElementX; i += step) {
		for (double j = startElementY; j < startElementY + endElementY; j += step) {
			temp.x = i;
			temp.y = j;
			temp.value = func(temp);
			if (min.value > temp.value) {
				min = temp;
			}
		}
	}
	return min;
}

int main(int argv, char** argc)
{
	int size, rank;
	if (MPI_Init(&argv, &argc) != MPI_SUCCESS)//Проверка на инициализацию
		return 1;
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS)//Получение размера коммуникатора
		return 2;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)//Получение текущего ранга 
		return 3;
	//if (size < 2)return 4;
	MPI_Status status;
	

	//Время
	double start = MPI_Wtime();



	//Работа с топологией
	MPI_Comm cart_comm;
	int cords[ndim];
	int dims[ndim] = { 0 };
	int periods[ndim] = { 0 };
	MPI_Dims_create(size, ndim, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 0, &cart_comm);
	MPI_Cart_coords(cart_comm, rank, ndim, cords);

	double step = hStart;

	//Выделение участков под каждый процесс
	int sizeX, remainsX;
	sizeX = ((right-left)/step )/dims[0];
	remainsX = (int)((right - left) / step) % dims[0];

	int* colX = new int[size];

	for (int i = 0; i < size; i++)colX[i] = sizeX;
	int counter = 0;
	
	while (remainsX > 0) {
		colX[counter++]++;
		if (counter >= size- 1)counter = 0;
		remainsX--;
	}

	//Нахождение стартовых элементов
	double startElementX = left;
	for (int i = 0; i < cords[0];i++) {
		startElementX += colX[i] * step;
	}
	point min;

	//Нахождение локального минимума
	min = findMin(startElementX, colX[rank] * step, bottom, top, step);

	//Нахождение глобального минимума
	double * buf = new double[3];
	cout << rank << " My value:" << min.value << endl;
	double* globalMinimum = new double[3];
	if (rank == 0) {
		double** results = new double*[size];
		results[0] = new double[3];
		results[0][0] = min.value;
		results[0][1] = min.x;
		results[0][2] = min.y;
		globalMinimum[0] = min.value;
		globalMinimum[1] = min.x;
		globalMinimum[2] = min.y;
		for (int i = 1; i < size; i++) {
			results[i] = new double[3];
			MPI_Recv(results[i], 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, cart_comm, &status);
		}
		
		for (int i = 1; i < size; i++) {
			if (results[i][0] < globalMinimum[0]) {
				globalMinimum[0] = results[i][0];
				globalMinimum[1] = results[i][1];
				globalMinimum[2] = results[i][2];
			}
		}
		for (int i = 1; i < size; i++) {
			MPI_Send(globalMinimum, 3, MPI_DOUBLE, i, 0, cart_comm);
		}
		
	}
	else {
		globalMinimum[0] = min.value;
		globalMinimum[1] = min.x;
		globalMinimum[2] = min.y;
		MPI_Send(globalMinimum, 3, MPI_DOUBLE, 0, 0, cart_comm);
		MPI_Recv(globalMinimum, 3, MPI_DOUBLE, 0, 0, cart_comm, &status);
	}


	step = hEnd;
	startElementX = globalMinimum[1]-hStart;
	for (int i = 0; i < cords[0]; i++) {
		startElementX += colX[i] * step;
	}

	//Нахождение локального минимума
	min = findMin(startElementX, colX[rank] * step, globalMinimum[2] - hStart, globalMinimum[2] + hStart, step);
	

	if (rank == 0) {
		
		double** results = new double* [size];
		results[0] = new double[3];
		results[0][0] = min.value;
		results[0][1] = min.x;
		results[0][2] = min.y;
		globalMinimum[0] = min.value;
		globalMinimum[1] = min.x;
		globalMinimum[2] = min.y;
		for (int i = 1; i < size; i++) {
			results[i] = new double[3];
			MPI_Recv(results[i], 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, cart_comm, &status);
		}

		for (int i = 1; i < size; i++) {
			if (results[i][0] < globalMinimum[0]) {
				globalMinimum[0] = results[i][0];
				globalMinimum[1] = results[i][1];
				globalMinimum[2] = results[i][2];
			}
		}
		
	}
	else {
		globalMinimum[0] = min.value;
		globalMinimum[1] = min.x;
		globalMinimum[2] = min.y;
		MPI_Send(globalMinimum, 3, MPI_DOUBLE, 0, 0, cart_comm);
	}

	if (rank == 0) {
		cout << "### x:" << globalMinimum[1] << " y:" << globalMinimum[2] << " value:" << globalMinimum[0] << endl;
		double end = MPI_Wtime();
		cout << "Time: " << end - start;
	}
	MPI_Finalize();
}