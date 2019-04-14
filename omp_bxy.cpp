#include <omp.h> 
#include<stdio.h>
#include <time.h>
static long num_steps = 100000;
double step;
#include<Windows.h>
#include<iostream>
#include<math.h>
using namespace std;

#define epx 0.5 * 1e-5
#define M 10000000 //最大迭代次数
#define N 5  //方程组的阶数
#define tread_count 4

int main() {
	clock_t start=clock();
	double a[N][N] = { { 0,0,0,0,0 },{ 0,8,7,0,0 },{ 0,6, 12, 5, 0 },{ 0,0, 4, 9, 3 },{ 0, 0, 0, 1, 2 } };
	double b[N] = { 0, 0, -2, 8,6 };
	double y[N]; //结果向量 
	double x[N];
	double norm;
	memset(y, 0, sizeof(y));
	int i, j;
	int flag;
	double result;

	for(int k=0;k<M;k++)
	{
		flag = 0;
# pragma omp parallel for
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
	clock_t finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
	//打印结果s
	for (i = 1; i < N; i++) {
		printf("%-4lf\t", y[i]);
	}
	printf("\nTime used:%lf s\n", duration);

	return 0;
}

