#include"cstdio.h"
#include"cmath.h"
#include"iostream.h"
#include"definer.h"
#include"functions.h"

int trade(double*b, int n, int m, double*d)
{
	if (m != 3 * n - 2)//首先判断是否满足三对角矩阵的结构条件
	{
		cout << "矩阵不满足三对角矩阵的条件！" << endl;
		return -2;
	}
	for (int k = 0; k <= n - 2; k++)
	{
		int j = k * 3;
		double s = b[j];
		if (fabs(s) + 1.0 == 1.0)//fabs()函数是求浮点数的绝对值，返回类型为double
		{
			cout << "分母为0、计算错误！" << endl;
			return -1;
		}
		b[j + 1] = b[j + 1] / s;//系数矩阵归一化
		d[k] = d[k] / s;//常数向量归一化
		b[j + 3] = b[j + 3] - b[j + 2] * b[j + 1];//系数矩阵消元
		d[k + 1] = d[k + 1] - b[j + 2] * d[k];
	}
	s = b[3 * n - 3];
	if (fabs(s) + 1.0 == 1.0)
	{
		cout << "分母为0，计算有误！" << endl;
		return -1;
	}
	d[n - 1] = d[n - 1] / s;//回带，解出x(n-1)
	for (k = n - 2; k >= 0; k--)
	{
		d[k] = d[k] - b[3 * k + 1] * d[k + 1];
	}
	return 2;

}