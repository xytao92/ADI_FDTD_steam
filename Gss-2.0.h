/*
调用GSS2.0求解稀疏矩阵
稀疏矩阵采用压缩列格式存储，下标从0开始

以列方向为标准方向。
变量定义：
ptr----->ptr的长度是N+1。ptr[i]记录第i列第一个非零元的位置,最后一个元素ptr[N]=nnz
ind----->每个非零元的行标
val------>每个非零元的数值


N=Nx或Ny或Nz
nnz=3*N-2

压缩列格式采用（int* ptr,int* ind,double *val）记录稀疏矩阵的非零元。
对于N阶有nnz个非零元的矩阵来说:
ind, val的长度是nnz。按列顺序记录每个非零元的行标和数值。
ptr的长度是N+1。ptr[i]记录第i列第一个非零元的位置,最后一个元素ptr[N]=nnz。
这样第i列的长度就是   ptr[i+1]-ptr[i]。

例如对于如下3阶的矩阵
1.0  0.0  5.0
0.0  3.0  6.0
2.0  4.0  7.0

假设压缩存储矩阵为B[],则有B[7]={1.0,2.0,3.0,4.0,5.0,6.0,7.0}按列优先进行存储

则 	ptr[4]={0,2,4,7};代表的是B[0],B[2],B[4],B[7]
ind[7]={0,2,1,2,0,1,2};
val[7]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

Copyright (c) 2005-present by YingShiChen.    All Rights Reserved.
Any problem,contact gsp@grusoft.com
*/
//#ifndef _GSS_2_0_
//#define _GSS_2_0_

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GSS_6_IMPORT __declspec(dllimport)		//申明下列函数从DLL库引入
#pragma comment( lib, "GSS_DLL_6D.lib")
#ifdef __cplusplus
extern "C" {
#endif   //让编译器采用C语言编译，防止的是C语言和C++语言的混合编译出现问题
//不用会出现这种错：错误	5	error LNK2019: 无法解析的外部符号 __imp__GSS_clear_ld，
//该符号在函数 "void __cdecl matel_gsscalc_bx(class Grid *,class Grid *)" (?matel_gsscalc_bx@@YAXPAVGrid@@0@Z) 中被引用...这是什么函数名，喂！

#define GRUS_MF_STATUS	0	
#define GRUS_OK	0

/*
	初始化

	int nRow,nCol	行数，列数
	int *ptr, int *ind, double *val	矩阵按压缩列格式存贮，(ptr,ind,val)为列指针，行下标和元素值	
	type	矩阵类型，定义如下。	
		0：		缺省值为0
		11：	对称正定矩阵。(ptr, ind, val)为矩阵下三角(含对角线)的数据。
		12：	对称不定矩阵。(ptr, ind, val)为矩阵下三角(含对角线)的数据。
	setting[32]		控制参数。

	返回值：
		初始化成功返回GRUS_OK,否则返回错误代码
*/
GSS_6_IMPORT	int GSS_init_ld(int nRow, int nCol, int* ptr, int* ind, double *val, int type, double *setting);

/*
	符号分解

	int nRow,nCol	行数，列数
	int *ptr, int *ind, double *val	矩阵按压缩列格式存贮，(ptr,ind,val)为列指针，行下标和元素值	

	返回值：
					分解成功返回求解器的指针,否则返回0x0
*/
GSS_6_IMPORT	void* GSS_symbol_ld(int nRow, int nCol, int* ptr, int* ind, double *val);

/*
	LU数值分解

	int nRow,nCol	行数，列数
	int *ptr, int *ind, double *val	矩阵按压缩列格式存贮，(ptr,ind,val)为列指针，行下标和元素值	
	void *hSolver		指向求解器的指针

	返回值：
					分解成功返回GRUS_OK, 否则返回错误代码
*/
GSS_6_IMPORT	int GSS_numeric_ld( int nRow,int nCol,int* ptr,int* ind,double *val,void *hSolver );

/*
	回代求解

	void *hSolver	指向求解器的指针
	int nRow,nCol	行数，列数
	int *ptr, int *ind, double *val	矩阵按压缩列格式存贮，(ptr,ind,val)为列指针，行下标和元素值	
	double *rhs		方程组的右端项，返回的解也存储在rhs中

	返回值：
					完成求解返回GRUS_OK, 否则返回错误代码
*/
GSS_6_IMPORT	int GSS_solve_ld( void *hSolver,int nRow,int nCol,int *ptr,int *ind,double *val,double *rhs );

/*
	释放求解器占用的内存

	void *hSolver	指向求解器的指针

	返回值：
					完成释放返回GRUS_OK, 否则返回错误代码
*/
GSS_6_IMPORT	int GSS_clear_ld( void* hSolver );

#ifdef __cplusplus
}
#endif
#pragma comment( lib, "GSS_DLL_6D.lib")
//#endif
