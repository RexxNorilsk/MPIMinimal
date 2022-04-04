﻿#include <iostream>
#include "mpi.h"
#include <Windows.h>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <fstream>
//Вариант 9
//mpiexec -n 4 "C:\Parallel 3.2\MinimalGrid\x64\Debug\MinimalGrid.exe"

//Топология
#define ndim 2

//Шаги
#define hCrash 3
#define hStart 0.01
#define hEnd 0.0001

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
	int dims[ndim] = { 0,0 };
	int periods[ndim] = { 0,0 };
	MPI_Dims_create(size, ndim, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 0, &cart_comm);
	MPI_Cart_coords(cart_comm, rank, ndim, cords);
	int sourceX, destX, sourceY, destY;
	MPI_Cart_shift(cart_comm, 0, 1, &sourceX, &destX);
	MPI_Cart_shift(cart_comm, 1, 1, &sourceY, &destY);

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
	int startElementX = left;
	for (int i = 0; i < cords[0];i++) {
		startElementX += colX[i];
	}
	int colY = top - bottom;
	int startElementY = bottom;
	point min;

		//Нахождение локального минимума
		min = findMin(startElementX, colX[rank] * step, startElementY, colY * step, step);

		//Нахождение глобального минимума
		double * buf = new double[3];
		cout << rank << " My value:" << min.value << endl;
		//Заменить на отправку на нулевой
		//MPI_Allreduce(&buf, &result, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
		if (rank == 0) {
			point* minList = new point[size];
			double** results = new double*[size];
			results[0] = new double[3];
			results[0][0] = min.value;
			results[0][1] = min.x;
			results[0][2] = min.y;
			for (int i = 1; i < size; i++) {
				results[i] = new double[3];
				MPI_Recv(results[i], 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, cart_comm, &status);
			}
			
			for (int i = 1; i < size; i++) {
				if (results[i][0] < min.value) {
					min.value = results[i][0];
					min.x = results[i][1];
					min.y = results[i][2];
				}
				if (results[i][0] == min.value) {
					point temp;
					temp.value = results[i][0];
					temp.x = results[i][1];
					temp.y = results[i][2];

					minList[0] = temp;
					minList[0] = temp;
					minList[0] = temp;
				}
			}
		}

		//Отправка данных
		double* data;
		if (result == buf) {
			data = new double[3];
			data[0] = min.x;
			data[1] = min.y;
			data[2] = min.value;
			for (int i = 0; i < size; i++) {
				if (i != rank) {
					MPI_Send(&data, 3, MPI_DOUBLE, i, 0, cart_comm);
					cout << rank << " To " << i << endl;
				}
			}
		}
		else {
			data = new double[3];
			MPI_Recv(&data, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, cart_comm, &status);
		}
		startElementX = data[0] - step;
		for (int i = 0; i < cords[0]; i++) {
			startElementX += colX[i];
		}
		startElementY = data[1] - step;
		step /= hCrash;

	if (rank == 0) {
		cout << "### x:" << min.x << " y:" << min.y << " value:" << min.value << endl;
		double end = MPI_Wtime();
		cout << "Time: " << end - start;
	}


	MPI_Finalize();
}