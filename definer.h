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

const double Z = 0.1; //单位m，波导长度，暂且设置的较长以便观察数据
const double X =  0.02286; //单位m，波导宽度，依据BJ100标准宽度
const double Y =  0.01016; //单位m, 波导高度，依据BJ100标准高度

const int Nz = 51; //长度方向网格数
const int Nx = 21;//宽度方向网格数
const int Ny = 11;//高度方向网格数

const double dz = Z/(Nz-1); //l方向空间步长
const double dx = X/(Nx-1); //w方向空间步长
const double dy = Y/(Ny-1); //h方向空间步长

const double dt = 1.0e-12; // 时间步长,  T=100dt ,T为输入源周期稳定性条件要求其小于2.3663e-12
const int STEPS = 600; //计算600个时间步长的数据

const double freq = 10.0e9; //入射场频率

const double T = 1/freq; //入射场周期10e-10s
const double pi = 3.14159265359;
const double c = 2.99792458e8;
const double omega = 2*pi*freq;//角速度//

const double epsl0 = 8.854e-12;//自由空间介电常数
const double mur0 = 4.0*pi*1.0e-7;//自由空间磁导率
const double k0 = omega*sqrt(mur0*epsl0);//自由空间波数//
const double kc = pi/X; //TE10模截止波数//
const double bate = sqrt(omega*omega*mur0*epsl0 - (pi/X)*(pi/X));
const double hm = 10.0;

const double epsl_z = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dx * dx)));
const double epsl_x = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dy * dy)));
const double epsl_y = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dz * dz)));

const double mur_z = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dx * dx)));
const double mur_x = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dy * dy)));
const double mur_y = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dz * dz)));
//定义输出文件
ofstream file_matle("result\\temp_matle_7.txt");
ofstream file_matle_plat("result\\platform_matle_7.txt");
ofstream file_mur1("result\\temp_mur1_7.txt");
ofstream file_mur1_plat("result\\platform_mur1_7.txt");

//class GSS
//{
//public:
//	int nRet = 0;//GSS函数的返回值
//	int N = Ny;
//	int nnz = 3 * N - 2;
//	int nRow = N;
//	int nCol = N;
//
//	int ptr[Ny + 1];
//	int ind[3 * Ny - 2];
//	double val[3 * Ny - 2];
//	double rhs[Ny];
//
//	void *hSolver = NULL;//求解器指针
//	double setting[32];
//	//for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
//	int type = 0;
//
//};

#endif