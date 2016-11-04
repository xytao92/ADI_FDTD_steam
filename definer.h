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

const int Nx = 501; //长度方向网格数
const int Ny = 201;//宽度方向网格数
const int Nz = 101;//高度方向网格数

const double  dx = 2e-4; //l方向空间步长
const double dy = 1.143e-4; //w方向空间步长
const double dz = 1.1016e-4; //h方向空间步长

const double dt = 1.0e-14; // 时间步长
const int STEPS = 12000; //计算6000个时间步长的数据

const double freq = 9e9; //入射场频率
const double pi = 3.1415926;
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



#endif