//writen by liuzhichao 
//2016/10/8
// i-1/2 ---> i-1;    i+1/2---> i;   i+1---> i+1半网格点的对应方式
#ifndef _DEFINER_
#define _DEFINER_


#include<cstdio>
#include <iostream>
#include <stdlib.h>
#include <fstream> 
#include <string>
#include <math.h>
#include <time.h>
#include<windows.h>

using namespace std;
//要求使用TE10模式作为强迫性激励源

const double X = 0.1; //单位m，波导长度，暂且设置的较长以便观察数据
const double Y = 0.02286; //单位m，波导宽度，依据BJ100标准宽度
const double Z = 0.01016; //单位m, 波导高度，依据BJ100标准高度

const int Nx = 51; //长度方向网格数
const int Ny = 21;//宽度方向网格数
const int Nz = 11;//高度方向网格数

const double  dx = 2e-4; //l方向空间步长
const double dy = 1.143e-4; //w方向空间步长
const double dz = 1.1016e-4; //h方向空间步长

const double dt = 1.0e-12; // 时间步长,  T=100dt ,T为输入源周期
const int STEPS = 600; //计算6000个时间步长的数据

const double freq = 1.0e10; //入射场频率
const double pi = 3.14159265359;
const double c = 2.99792458e8;
const double omega = 2 * pi*freq;//角速度//

const double epsl0 = 8.854e-12;//自由空间介电常数
const double mur0 = 4.0*pi*1.0e-7;//自由空间导磁率
const double  k0 = omega*sqrt(mur0*epsl0);//自由空间波数//
const double kc = pi / Y; //TE10模截止波数//


const double epsl_x = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dy * dy)));
const double epsl_y = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dz * dz)));
const double epsl_z = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dx * dx)));

const double mur_x = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dy * dy)));
const double mur_y = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dz * dz)));
const double mur_z = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dx * dx)));


class GSS
{
public:
	int nRet = 0;//GSS函数的返回值
	int N = Ny;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Ny + 1];
	int ind[3 * Ny - 2];
	double val[3 * Ny - 2];
	double rhs[Ny];

	void *hSolver = NULL;//求解器指针
	double setting[32];
	//for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
	int type = 0;

};

#endif