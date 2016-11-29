/*
Writen By LiuZhiChao 
on time: 2016.11.22
All Rights Receved
*/

#ifndef _GSS_ADI_
#define _GSS_ADI_

#include "definer.h"
#include "grid.h"
#include"functions.h"
#include"Gss-2.0.h"
//-----------------------------------------------------------------------函数声明
void matel_gsscalc_ex(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void matel_gsscalc_ey(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void matel_gsscalc_ez(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void matel_gsscalc_bx(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void matel_gsscalc_by(Grid* halfgrid_before, Grid* halfgrid_now, int step);
void matel_gsscalc_bz(Grid* halfgrid_before, Grid* halfgrid_now, int step);

void adi_fdtd_leapforg_matel_GSS(Grid* halfgrid_before, Grid* halfgrid_now)

{
	//由于是理想金属的，那么可以在边界处赋值为0，可以直接用gss求解
	system("mkdir result");
	ofstream file("result\\GSS_matel_steam2_1.txt");//用于保存结果

	//*******计算TE10模******//注意边界条件的问题，不处理周围的四个面
	//PART1---- 计算电场//

	int step = 0;//计算时间步长
	while (step < STEPS)
	{

		inject_field(halfgrid_before,halfgrid_now, step);
		//---------------------------------------------------------------------------------------计算电场
		//基本思想：分别计算六个参量全网格的值，由于计算不需要当前时刻的值，可以分步全局计算

		matel_gsscalc_ex(halfgrid_before, halfgrid_now,step);
		matel_gsscalc_ey(halfgrid_before, halfgrid_now,step);
		matel_gsscalc_ez(halfgrid_before, halfgrid_now,step);
		
		//---------------------------------------------------------------------------------------计算磁场
		matel_gsscalc_bx(halfgrid_before, halfgrid_now,step);
		matel_gsscalc_by(halfgrid_before, halfgrid_now,step);
	    matel_gsscalc_bz(halfgrid_before, halfgrid_now,step);

		int result_x = 0;
		int result_y = 10;
		int result_z = 5;
		
		file << step << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ex << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ey << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ez << '\t';
		file << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bx << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].by << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bz << '\t';
		file << '\n';
		cout << "Step--- " << step << " ---has finished." << endl;
		step++;
	}//while
}//函数结尾

void matel_gsscalc_ex(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	double aex = (-1 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dy)*(1 / dy));
	double bex = (epsl0 + 2 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dy)*(1 / dy));
	double cex = aex;

	int nRet;//GSS函数的返回值
	int N = Ny-1;
	
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Ny];
	int ind[3 * (Ny-1) - 2];
	double val[3 * (Ny-1) - 2];
	double rhs[Ny-1];

	//----------------------------------------------------------------------------test
	/*int ptr1[4] = { 0, 2, 5, 7 };
	int ind1[7] = {0,1,0,1,2,1,2};
	double val1[7] = {bex,cex,aex,bex,cex,aex,bex};

	double temp1 = bex*halfgrid_before[0].ex + cex*halfgrid_before[Nz].ex
					+ dt*((halfgrid_before[0].bz) / dy - (halfgrid_before[0].by) / dz);
	double temp2 = bex*halfgrid_before[Nz].ex + cex*halfgrid_before[2*Nz].ex
				+ dt*((halfgrid_before[Nz].bz) / dy - (halfgrid_before[Nz].by) / dz);
	double temp3 = bex*halfgrid_before[2*Nz].ex + cex*halfgrid_before[3*Nz].ex
					+ dt*((halfgrid_before[2*Nz].bz) / dy - (halfgrid_before[2*Nz].by) / dz);
		
	double rhs1[3] = {temp1,temp2,temp3};*/
	//----------------------------------------------------------------------------test

	void *hSolver=NULL;//求解器指针
	double setting[32];
	for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
	int type = 0;
	//----------------------------------------------------------------------------
	//---------------------------------------------------------------处理ptr数组
	ptr[0] = 0;
	ptr[1] = 2;
	ptr[N] = 3 * N - 2;
	for (int i = 2; i < N; i++)
	{
		ptr[i] = ptr[i - 1] + 3;
	}
	
	//--------------------------------------------------------------处理ind数组
	ind[0] = 0;
	ind[1] = 1;
	ind[3 * N - 3] = N - 1;
	ind[3 * N - 4] = N - 2;
	for (int i = 2,  j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//--------------------------------------------------------------val数组处理
	val[0] = bex;
	val[2] = cex;
	val[3 * N - 3] = bex;
	val[3 * N - 5] = aex;
	for (int i = 1,  j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = aex;
		val[i + 2] = bex;
		val[i + 4] = cex;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;
    //------------------------------------------------------------生成系数矩阵

	for (int i = 0; i < Nx-1; i++)
	{
		for (int k = 0; k < Nz-1; k++)
		{
			for (int j = 0; j < Ny-1; j++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (j == 0 && k != 0)
					rhs[j] = bex*halfgrid_before[i*Ny*Nz + j*Nz + k].ex + cex*halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bz) / dy - (halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].by) / dz);
				else if (k == 0 && j != 0)
					rhs[j] = aex*halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex + bex*halfgrid_before[i*Ny*Nz + j*Nz + k].ex + cex*halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bz) / dy - (halfgrid_before[i*Ny*Nz + j*Nz + k].by) / dz);
				else if (j == 0 && k == 0)
					rhs[j] = bex*halfgrid_before[i*Ny*Nz + j*Nz + k].ex + cex*halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bz) / dy - (halfgrid_before[i*Ny*Nz + j*Nz + k].by) / dz);
				else if (j!=0&&k!=0)
					rhs[j] = aex*halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex + bex*halfgrid_before[i*Ny*Nz + j*Nz + k].ex + cex*halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bz) / dy - (halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].by) / dz);
			}
			
		
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK)	{
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL)	{
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK)	{
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);
			
			for (int j1 = 0; j1 < Ny-1; j1++)
			{//保存结果
				halfgrid_now[i*Ny*Nz + j1*Nz + k].ex = rhs[j1];
				halfgrid_before[i*Ny*Nz + j1*Nz + k].ex = halfgrid_now[i*Ny*Nz + j1*Nz + k].ex;
			}
			
			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
		
	}
	
}

void matel_gsscalc_ey(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{

	double aey = (-1 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dz)*(1 / dz));
	double bey = (epsl0 + 2 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dz)*(1 / dz));
	double cey = aey;

	int nRet;//GSS函数的返回值
	int N = Nz-1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nz];
	int ind[3 * (Nz-1) - 2];
	double val[3 * (Nz-1) - 2];
	double rhs[Nz-1];

	void *hSolver = NULL;//求解器指针
	double setting[32];
	for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
	int type = 0;

	//处理ptr数组
	ptr[0] = 0;
	ptr[1] = 2;
	ptr[N] = 3 * N - 2;
	for (int i = 2; i < N; i++)
	{
		ptr[i] = ptr[i - 1] + 3;
	}
	//处理ind数组
	ind[0] = 0;
	ind[1] = 1;
	ind[3 * N - 3] = N - 1;
	ind[3 * N - 4] = N - 2;
	for (int i = 2,  j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//val数组处理
	val[0] = bey;
	val[2] = cey;
	val[3 * N - 3] = bey;
	val[3 * N - 5] = aey;
	for (int i = 1,  j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = aey;
		val[i + 2] = bey;
		val[i + 4] = cey;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;
	
	for (int i = 0; i < Nx-1; i++)
	{
		for (int j = 0; j < Ny-1; j++)
		{
			for (int k = 0; k < Nz-1; k++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (k == 0 && i != 0)
					rhs[k] = bey*halfgrid_before[i*Ny*Nz + j*Nz + k].ey + cey*halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bx) / dz - (halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].bz) / dx);
				else if (i == 0 && k != 0)
					rhs[k] = aey*halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey + bey*halfgrid_before[i*Ny*Nz + j*Nz + k].ey + cey*halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].bx) / dz - (halfgrid_before[i*Ny*Nz + j*Nz + k].bz) / dx);
				else if (k == 0 && i == 0)
					rhs[k] = bey*halfgrid_before[i*Ny*Nz + j*Nz + k].ey + cey*halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bx) / dz - (halfgrid_before[i*Ny*Nz + j*Nz + k].bz) / dx);
				else if (k != 0 && i != 0)
					rhs[k] = aey*halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey + bey*halfgrid_before[i*Ny*Nz + j*Nz + k].ey + cey*halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].bx) / dz - (halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].bz) / dx);
			}//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK)	{
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL)	{
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK)	{
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int k1 = 0; k1 < Nz-1; k1++){//保存结果
				halfgrid_now[i*Ny*Nz + j*Nz + k1].ey = rhs[k1];
				halfgrid_before[i*Ny*Nz + j*Nz + k1].ey = halfgrid_now[i*Ny*Nz + j*Nz + k1].ey;
			}

			if (hSolver != NULL)	
				GSS_clear_ld(hSolver);
		}
	}

}

void matel_gsscalc_ez(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	double aez = (-1 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dx)*(1 / dx));
	double bez = (epsl0 + 2 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dx)*(1 / dx));
	double cez = aez;

	int nRet;//GSS函数的返回值
	int N = Nx-1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nx];
	int ind[3 * (Nx-1) - 2];
	double val[3 * (Nx-1) - 2];
	double rhs[Nx-1];

	void *hSolver = NULL;//求解器指针
	double setting[32];
	for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
	int type = 0;

	//处理ptr数组
	ptr[0] = 0;
	ptr[1] = 2;
	ptr[N] = 3 * N - 2;
	for (int i = 2; i < N; i++)
	{
		ptr[i] = ptr[i - 1] + 3;
	}
	//处理ind数组
	ind[0] = 0;
	ind[1] = 1;
	ind[3 * N - 3] = N - 1;
	ind[3 * N - 4] = N - 2;
	for (int i = 2,  j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//val数组处理
	val[0] = bez;
	val[2] = cez;
	val[3 * N - 3] = bez;
	val[3 * N - 5] = aez;
	for (int i = 1,  j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = aez;
		val[i + 2] = bez;
		val[i + 4] = cez;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;
	
	for (int j = 0; j < Ny-1; j++)
	{
		for (int k = 0; k < Nz-1; k++)
		{
			for (int i = 0; i < Nx-1; i++)
			{
				////if ( i == 0&& j == 10&&k == 5)//加点源位置，保证不会被覆盖
				////	halfgrid_before[i*Ny*Nz + j*Nz + k].ez = 100 * sin(omega*step*dt);//100为设置值//
				//首先处理矩阵方程的右端项rhs数组
			     if (i == 0 && j != 0)
					rhs[i] =  bez*halfgrid_before[i*Ny*Nz + j*Nz + k].ez + cez*halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].by) / dx - (halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bx) / dy);
				else if (j == 0 && i != 0)
					rhs[i] = aez*halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez + bez*halfgrid_before[i*Ny*Nz + j*Nz + k].ez + cez*halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].by) / dx - (halfgrid_before[i*Ny*Nz + j*Nz + k].bx) / dy);
				else if (i == 0 && j == 0)
					rhs[i] =bez*halfgrid_before[i*Ny*Nz + j*Nz + k].ez + cez*halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].by) / dx - (halfgrid_before[i*Ny*Nz + j*Nz + k].bx) / dy);
				else if (i != 0 && j != 0)
					rhs[i] = aez*halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez + bez*halfgrid_before[i*Ny*Nz + j*Nz + k].ez + cez*halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez
					+ dt*((halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].by) / dx - (halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bx) / dy);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK)	{
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL)	{
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK)	{
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int i1 = 0; i1 < Nx-1; i1++)
			{//保存结果
				if (i1 == 0)
				{
					//halfgrid_before[j*Nz + k].ez = 100 * sin((pi / Y)*j*dy)*sin(omega*step*dt);//100为设置值//115
					halfgrid_now[j*Nz + k].ez = 100 * sin((pi / Y)*j*dy)*sin(omega*step*dt);//100为设置值//115
				}
				else
				{
					halfgrid_now[i1*Ny*Nz + j*Nz + k].ez = rhs[i1];

				     halfgrid_before[i1*Ny*Nz + j*Nz + k].ez = halfgrid_now[i1*Ny*Nz + j*Nz + k].ez;
				}
				
			}
					
			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
	
}


void matel_gsscalc_bx(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	double ahx = (-1 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dy)*(1 / dy));
	double bhx = (mur0 + 2 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dy)*(1 / dy));
	double chx = ahx;

	int nRet;//GSS函数的返回值
	int N = Ny-1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Ny];
	int ind[3 * (Ny-1) - 2];
	double val[3 * (Ny-1) - 2];
	double rhs[Ny-1];

	void *hSolver = NULL;//求解器指针
	double setting[32];
	for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
	int type = 0;

	//处理ptr数组
	ptr[0] = 0;
	ptr[1] = 2;
	ptr[N] = 3 * N - 2;
	for (int i = 2; i < N; i++)
	{
		ptr[i] = ptr[i - 1] + 3;
	}
	//处理ind数组
	ind[0] = 0;
	ind[1] = 1;
	ind[3 * N - 3] = N - 1;
	ind[3 * N - 4] = N - 2;
	for (int i = 2,  j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//val数组处理
	val[0] = bhx;
	val[2] = chx;
	val[3 * N - 3] = bhx;
	val[3 * N - 5] = ahx;
	for (int i = 1,  j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = ahx;
		val[i + 2] = bhx;
		val[i + 4] = chx;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;
	
	for (int i = 0; i < Nx-1; i++)
	{
		for (int k = 0; k < Nz-1; k++)
		{
			for (int j = 0; j < Ny-1; j++)
			{
				//首先处理矩阵方程的右端项rhs数组
 				if (j == 0)
					rhs[j] = bhx*halfgrid_before[i*Ny*Nz + j*Nz + k].bx + chx*halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].bx
					+ dt*((halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ey - halfgrid_now[i*Ny*Nz + j*Nz + k].ey) / dz - (halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ez - halfgrid_now[i*Ny*Nz + j*Nz + k].ez) / dy);
				else
					rhs[j] = ahx*halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bx + bhx*halfgrid_before[i*Ny*Nz + j*Nz + k].bx + chx*halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].bx
					+ dt*((halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ey - halfgrid_now[i*Ny*Nz + j*Nz + k].ey) / dz - (halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ez - halfgrid_now[i*Ny*Nz + j*Nz + k].ez) / dy);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK)	{
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL)	{
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK)	{
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int j1 = 0; j1 < Ny-1; j1++){//保存结果
				halfgrid_now[i*Ny*Nz + j1*Nz + k].bx = rhs[j1];
				halfgrid_before[i*Ny*Nz + j1*Nz + k].bx = halfgrid_now[i*Ny*Nz + j1*Nz + k].bx;
			}

			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
	
}

void matel_gsscalc_by(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{

	double ahy = (-1 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dz)*(1 / dz));
	double bhy = (mur0 + 2 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dz)*(1 / dz));
	double chy = ahy;

	int nRet;//GSS函数的返回值
	int N = Nz-1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nz];
	int ind[3 * (Nz-1) - 2];
	double val[3 * (Nz-1) - 2];
	double rhs[Nz-1];

	void *hSolver = NULL;//求解器指针
	double setting[32];
	for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
	int type = 0;

	//处理ptr数组
	ptr[0] = 0;
	ptr[1] = 2;
	ptr[N] = 3 * N - 2;
	for (int i = 2; i < N; i++)
	{
		ptr[i] = ptr[i - 1] + 3;
	}
	//处理ind数组
	ind[0] = 0;
	ind[1] = 1;
	ind[3 * N - 3] = N - 1;
	ind[3 * N - 4] = N - 2;
	for (int i = 2,  j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//val数组处理
	val[0] = bhy;
	val[2] = chy;
	val[3 * N - 3] = bhy;
	val[3 * N - 5] = ahy;
	for (int i = 1,  j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = ahy;
		val[i + 2] = bhy;
		val[i + 4] = chy;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;
	
	for (int i = 0; i < Nx-1; i++)
	{
		for (int j = 0; j < Ny-1; j++)
		{
			for (int k = 0; k < Nz-1; k++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (k == 0)
					rhs[k] = bhy*halfgrid_before[i*Ny*Nz + j*Nz + k].by + chy*halfgrid_before[i*Ny*Nz + j*Nz + k + 1].by
					+ dt*((halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ez - halfgrid_now[i*Ny*Nz + j*Nz + k].ez) / dx - (halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ex - halfgrid_now[i*Ny*Nz + j*Nz + k].ex) / dz);
				else
					rhs[k] = ahy*halfgrid_before[i*Ny*Nz + j*Nz + k - 1].by + bhy*halfgrid_before[i*Ny*Nz + j*Nz + k].by + chy*halfgrid_before[i*Ny*Nz + j*Nz + k + 1].by
					+ dt*((halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ez - halfgrid_now[i*Ny*Nz + j*Nz + k].ez) / dx - (halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ex - halfgrid_now[i*Ny*Nz + j*Nz + k].ex) / dz);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK)	{
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL)	{
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK)	{
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int k1 = 0; k1 < Nz-1; k1++){//保存结果
				halfgrid_now[i*Ny*Nz + j*Nz + k1].by = rhs[k1];
				halfgrid_before[i*Ny*Nz + j*Nz + k1].by = halfgrid_now[i*Ny*Nz + j*Nz + k1].by;
			}

			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
	
}

void matel_gsscalc_bz(Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	double ahz = (-1 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dx)*(1 / dx));
	double bhz = (mur0 + 2 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dx)*(1 / dx));
	double chz = ahz;

	int nRet;//GSS函数的返回值
	int N = Nx-1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nx];
	int ind[3 * (Nx-1) - 2];
	double val[3 * (Nx-1) - 2];
	double rhs[Nx-1];

	void *hSolver = NULL;//求解器指针
	double setting[32];
	for (int i = 0; i < 32; i++)	setting[i] = 0.0;//配置参数初始化
	int type = 0;

	//处理ptr数组
	ptr[0] = 0;
	ptr[1] = 2;
	ptr[N] = 3 * N - 2;
	for (int i = 2; i < N; i++)
	{
		ptr[i] = ptr[i - 1] + 3;
	}
	//处理ind数组
	ind[0] = 0;
	ind[1] = 1;
	ind[3 * N - 3] = N - 1;
	ind[3 * N - 4] = N - 2;
	for (int i = 2,  j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//val数组处理
	val[0] = bhz;
	val[2] = chz;
	val[3 * N - 3] = bhz;
	val[3 * N - 5] = ahz;
	for (int i = 1,  j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = ahz;
		val[i + 2] = bhz;
		val[i + 4] = chz;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;
	
	for (int j = 0; j < Ny-1; j++)
	{
		for (int k = 0; k < Nz-1; k++)
		{
			for (int i = 0; i < Nx-1; i++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (i == 0)
					rhs[i] = bhz*halfgrid_before[i*Ny*Nz + j*Nz + k].bz + chz*halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].bz
					+ dt*((halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ex - halfgrid_now[i*Ny*Nz + j*Nz + k].ex) / dy - (halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ey - halfgrid_now[i*Ny*Nz + j*Nz + k].ey) / dx);
				else
					rhs[i] = ahz*halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].bz + bhz*halfgrid_before[i*Ny*Nz + j*Nz + k].bz + chz*halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].bz
					+ dt*((halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ex - halfgrid_now[i*Ny*Nz + j*Nz + k].ex) / dy - (halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ey - halfgrid_now[i*Ny*Nz + j*Nz + k].ey) / dx);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK)	{
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL)	{
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK)	{
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int i1 = 0; i1 < Nx-1; i1++){//保存结果
				halfgrid_now[i1*Ny*Nz + j*Nz + k].bz = rhs[i1];
				halfgrid_before[i1*Ny*Nz + j*Nz + k].bz = halfgrid_now[i1*Ny*Nz + j*Nz + k].bz;
			}

			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
	
}

#endif