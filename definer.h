//writen by liuzhichao 
//2016/10/8
// i-1/2 ---> i-1;    i+1/2---> i;   i+1---> i+1半网格点的对应方式

//2017.03.29检查过参数，没有问题 by 刘智超
//2017.03.30问过金老师，修改网格的尺寸与dz一致,实验发现，结果不错。经验是不能让网格的大小三边差别太大。

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
/*不一定非要用标准的BJ100波导去计算，可取整数长度
*实验一：Z=0.4,X=0.02,Y=0.01;Nz=401,Nx=21,Ny=11;结果，波形稳定，但是感觉在逐步达到理论幅值，不知后续是不是会一直增加下去。
*实验二：Z=0.4,X=0.02,Y=0.01;Nz=201,Nx=21,Ny=11;目的看是否dz要是dz,dx的倍数关系。结果：波形稳定，但是感觉在逐步达到理论幅值，不知后续是不是会一直增加下去。
*实验证明：只要保证dz不变，等比例延长计算结构不改变计算结果。
*实验三：测试另外的点，看是否依然符合的比较好。结果，是的。
*实验四：1011*101*51，测试区域大小不变，效果不好，等比例增大网格数是不行的
*实验五：211*21*11，测试区域Z方向不变，另外两个都缩小十倍，400步长看不出效果，800计算条件下，误差更大，而且不均匀

*/
//04.07----网格的尺度调节到跟dz一样大小，相当于网格加密了100倍，测试结果可以

const double Z = 0.4*1.05; //单位m，波导长度，增加十个网格长度以增加十层CPML,要保持dz不变则在增加CPML的十层时也要相应增加整个结构的长度
const double X = 0.2; //单位m，波导宽度，依据BJ100标准宽度
const double Y = 0.1; //单位m, 波导高度，依据BJ100标准高度

const int Nz = 211; //长度方向网格数 201，再加上吸收边界的10层，为211层
const int Nx = 101;//宽度方向网格数
const int Ny = 51;//高度方向网格数
const int s0 = 199;//开始加吸收边界的位置
const int Nz_size = Nz - s0 - 2;


const double dz = Z / (Nz - 1); //z方向网格大小
const double dx = X / (Nx - 1); //x方向网格大小
const double dy = Y / (Ny - 1); //y方向网格大小

const double dt = 0.5*3.85167e-12; // CFLN=0.5,时间步长
const int STEPS = 400; //计算400个时间步长的数据

const double freq = 10.0e9; //入射场频率

const double T = 1/freq; //入射场周期1.0e-10s
const double pi = 3.14159265359;
const double c = 2.99792458e8;
const double omega = 2*pi*freq;//角速度//

const double e = 2.71828182845;
const double epsl0 = 8.854187817e-12;//自由空间介电常数
const double mur0 = 1.2566370614e-6;//自由空间磁导率
//const double k0 = omega*sqrt(mur0*epsl0);//自由空间波数//
const double kc = pi/X; //TE10模截止波数//
const double bate = sqrt(omega*omega*mur0*epsl0 - (pi/X)*(pi/X));
const double hm = 10.0;

//const double epsl_z = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dx * dx)));
//const double epsl_x = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dy * dy)));
//const double epsl_y = (epsl0 - (dt / 2) * (dt / 2) * ((1 / mur0) / (dz * dz)));
//
//const double mur_z = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dx * dx)));
//const double mur_x = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dy * dy)));
//const double mur_y = (mur0 - (dt / 2) * (dt / 2) * ((1 / epsl0) / (dz * dz)));
//定义输出文件
const string matle_filepath = "result\\matle_400_N0.5.txt";
const string matle_p_filepath = "result\\platform_matle_1500_N0.5.txt";

const string cpml_filepath = "result\\cpml_400_N0.5.txt";
const string cpml_p_filepath = "result\\platform_cpml_400_N0.5.txt";

const string theor_val_filepath = "result\\theor_val\\source_400_N0.5.txt";

#endif