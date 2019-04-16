#include <omp.h> 
#include<stdio.h>
#include <time.h>
static long num_steps = 100000;
double step;
#include<Windows.h>
#include<iostream>
#include<math.h>
using namespace std;

#define epx 0.01 * 1e-8
#define M 1000 //最大迭代次数
#define N 5  //方程组的阶数
#define tread_count 4

double a[N][N] = { { 0,8,5,2,1 },{ 0,8,7,0,0 },{ 0,6, 12, 5, 0 },{ 0,0, 4, 9, 3 },{ 3, 2, 0, 1, 2 } };
double b[N] = { 0, 0, -2, 8,6 };
double y[N]; //结果向量 
double x[N];
double norm;
int i, j;
int flag;
double result;

double para() {

	double start = omp_get_wtime();
	for (int k = 0; k<M; k++)
	{
		flag = 0;
#pragma omp parallel for
		for (i = 0; i < N; i++) {
			x[i] = y[i];
		}

#pragma omp parallel for
		for (i = 1; i < N; i++) {
			result = 0;
			for (j = 1; j < N; j++) {
				if (i != j) {
					result += a[i][j] * y[j];
				}
			}
			y[i] = (b[i] - result) / a[i][i]; //求出解空间
		}

		//求y -x 的范数
		norm = 0.0;
#pragma omp parallel for
		for (i = 1; i < N; i++) {
			if (fabs(x[i] - y[i]) > norm) {
				norm = fabs(x[i] - y[i]);
			}
		}

		if (norm < epx) {
			break;
		}
	}
	double finish = omp_get_wtime();
	double d = (finish - start);
	return d;

}

double serial() {

	double start = omp_get_wtime();
	for (int k = 0; k<M; k++)
	{
		flag = 0;
		for (i = 0; i < N; i++) {
			x[i] = y[i];
		}

		for (i = 1; i < N; i++) {
			result = 0;
			for (j = 1; j < N; j++) {
				if (i != j) {
					result += a[i][j] * y[j];
				}
			}
			y[i] = (b[i] - result) / a[i][i]; //求出解空间
		}

		//求y -x 的范数
		norm = 0.0;
		for (i = 1; i < N; i++) {
			if (fabs(x[i] - y[i]) > norm) {
				norm = fabs(x[i] - y[i]);
			}
		}

		if (norm < epx) {
			break;
		}
	}
	double finish = omp_get_wtime();
	double d = (finish - start);

	return d;

}


int main() {
	memset(y, 0, sizeof(y));
	//打印结果s
	double ser = serial();
	double par = para();
	for (i = 1; i < N; i++) {
		printf("%-4lf\t", y[i]);
	}
	printf("\n Speed Up:%lf s\n", ser/par);

	return 0;
}

