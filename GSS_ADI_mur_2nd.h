
/*
Writen By LiuZhiChao
on time: 2016.11.22
All Rights Receved
*/

#ifndef _GSS_ADI_MUR_2ND
#define _GSS_ADI_MUR_2ND

#include "definer.h"
#include "grid.h"
#include"functions.h"
#include"Gss-2.0.h"
//-----------------------------------------------------------------------函数声明
void mur2_gsscalc_ez(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step);
void mur2_gsscalc_ex(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step);
void mur2_gsscalc_ey(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step);

void mur2_gsscalc_bz(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step);
void mur2_gsscalc_bx(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step);
void mur2_gsscalc_by(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step);  

void adi_fdtd_leapforg_mur2_GSS(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now)

{
	//由于是理想金属的，那么可以在边界处赋值为0，可以直接用gss求解
	//system("mkdir result");
	ofstream file("result\\temp_mur_2nd.txt");//用于保存结果
	ofstream file2("result\\platform_ex_mur_2nd.txt");
	ofstream file3("result\\platform_ey_mur_2nd.txt");
	ofstream file4("result\\platform_ez_mur_2nd.txt");

	ofstream file5("result\\platform_bx_mur_2nd.txt");
	ofstream file6("result\\platform_by_mur_2nd.txt");
	ofstream file7("result\\platform_bz_mur_2nd.txt");
	//*******计算TE10模******//-------------------------注意边界条件的问题，不处理周围的四个面
	//PART1---- 计算电场//

	int step = 0;//计算时间步长
	while (step < STEPS)
	{

		inject_field(halfgrid_before, halfgrid_now, step);
		//---------------------------------------------------------------------------------------计算电场
		//基本思想：分别计算六个参量全网格的值，由于计算不需要当前时刻的值，可以分步全局计算

		mur2_gsscalc_ez(halfgrid_beforeX2,halfgrid_before, halfgrid_now, step);
		mur2_gsscalc_ex(halfgrid_beforeX2,halfgrid_before, halfgrid_now, step);
		mur2_gsscalc_ey(halfgrid_beforeX2,halfgrid_before, halfgrid_now, step);

		//---------------------------------------------------------------------------------------计算磁场
		mur2_gsscalc_bz(halfgrid_beforeX2,halfgrid_before, halfgrid_now, step);
		mur2_gsscalc_bx(halfgrid_beforeX2,halfgrid_before, halfgrid_now, step);
		mur2_gsscalc_by(halfgrid_beforeX2,halfgrid_before, halfgrid_now, step);

		int result_z = 5;
		int result_x = 10;
		int result_y = 5;//1270

		file << step << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ex << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ey << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ez << '\t';
		file << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].bx << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].by << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].bz << '\t';
		file << '\n';

		cout << "Step--- " << step << " ---has finished." << endl;
		step++;
	}//while
	for (int i = 0; i < Nx; i++)//输出一个横截面的数据
	{
		for (int k = 0; k < Nz; k++)
		{
			file2 << k << '\t' << i << '\t' << halfgrid_now[k*Nx*Ny + i*Ny + 5].ex << endl;
			file3 << k << '\t' << i << '\t' << halfgrid_now[k*Nx*Ny + i*Ny + 5].ey << endl;
			file4 << k << '\t' << i << '\t' << halfgrid_now[k*Nx*Ny + i*Ny + 5].ez << endl;

			file5 << k << '\t' << i << '\t' << halfgrid_now[k*Nx*Ny + i*Ny + 5].bx << endl;
			file6 << k << '\t' << i << '\t' << halfgrid_now[k*Nx*Ny + i*Ny + 5].by << endl;
			file7 << k << '\t' << i << '\t' << halfgrid_now[k*Nx*Ny + i*Ny + 5].bz << endl;
		}

	}

}//函数结尾

void mur2_gsscalc_ez(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	//cout << "源的Ez值---------------->" << halfgrid_now[10 * Ny + 5].ez << endl;
	double aez = (-1 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dx)*(1 / dx));
	double bez = (epsl0 + 2 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dx)*(1 / dx));
	double cez = aez;

	int nRet;//GSS函数的返回值
	int N = Nx - 1;

	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nx];
	int ind[3 * (Nx - 1) - 2];
	double val[3 * (Nx - 1) - 2];
	double rhs[Nx - 1];

	void *hSolver = NULL;//求解器指针
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
	for (int i = 2, j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//--------------------------------------------------------------val数组处理
	val[0] = bez;
	val[2] = cez;
	val[3 * N - 3] = bez;
	val[3 * N - 5] = aez;
	for (int i = 1, j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = aez;
		val[i + 2] = bez;
		val[i + 4] = cez;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;
	//------------------------------------------------------------生成系数矩阵

	for (int k = 1; k < Nz - 1; k++)
	{
		for (int j = 0; j < Ny - 1; j++)
		{
			for (int i = 0; i < Nx - 1; i++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (i == 0 && j != 0)
					rhs[i] = bez*halfgrid_before[k*Nx*Ny + i*Ny + j].ez + cez*halfgrid_before[k*Nx*Ny + (i + 1)*Ny + j].ez
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].by) / dx - (halfgrid_before[k*Nx*Ny + i*Ny + j].bx - halfgrid_before[k*Nx*Ny + i*Ny + j - 1].bx) / dy);
				else if (j == 0 && i != 0)
					rhs[i] = aez*halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j].ez + bez*halfgrid_before[k*Nx*Ny + i*Ny + j].ez + cez*halfgrid_before[k*Nx*Ny + (i + 1)*Ny + j].ez
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].by - halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j].by) / dx - (halfgrid_before[k*Nx*Ny + i*Ny + j].bx) / dy);
				else if (i == 0 && j == 0)
					rhs[i] = bez*halfgrid_before[k*Nx*Ny + i*Ny + j].ez + cez*halfgrid_before[k*Nx*Ny + (i + 1)*Ny + j].ez
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].by) / dx - (halfgrid_before[k*Nx*Ny + i*Ny + j].bx) / dy);
				else if (i != 0 && j != 0)
					rhs[i] = aez*halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j].ez + bez*halfgrid_before[k*Nx*Ny + i*Ny + j].ez + cez*halfgrid_before[k*Nx*Ny + (i + 1)*Ny + j].ez
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].by - halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j].by) / dx - (halfgrid_before[k*Nx*Ny + i*Ny + j].bx - halfgrid_before[k*Nx*Ny + i*Ny + j - 1].bx) / dy);
			}


			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK) {
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL) {
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK) {
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int i1 = 0; i1 < Nx - 1; i1++)
			{
				

				if (j == 0 || j == Ny - 2 || i1 == 0 || i1 == Nx - 2)//处理ex的边界问题，在四个面的位置应该为0,所有平面都进行处理
				{
					rhs[i1] = 0.0;
				}

				halfgrid_beforeX2[k*Nx*Ny + i1*Ny + j].ez = halfgrid_before[k*Nx*Ny + i1*Ny + j].ez;
				halfgrid_before[k*Nx*Ny + i1*Ny + j].ez = halfgrid_now[k*Nx*Ny + i1*Ny + j].ez;
				halfgrid_now[k*Nx*Ny + i1*Ny + j].ez = rhs[i1];				
			}
			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
}

void mur2_gsscalc_ex(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	//cout << "源的Ex值---------------->" << halfgrid_now[10 * Ny + 5].ex << endl;
	double aex = (-1 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dy)*(1 / dy));
	double bex = (epsl0 + 2 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dy)*(1 / dy));
	double cex = aex;

	int nRet;//GSS函数的返回值
	int N = Ny - 1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Ny];
	int ind[3 * (Ny - 1) - 2];
	double val[3 * (Ny - 1) - 2];
	double rhs[Ny - 1];

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
	for (int i = 2, j = 0; i + 2 < 3 * N - 4; i = i + 3)
	{
		ind[i] = j;
		ind[i + 1] = j + 1;
		ind[i + 2] = j + 2;
		j++;
	}
	//val数组处理
	val[0] = bex;
	val[2] = cex;
	val[3 * N - 3] = bex;
	val[3 * N - 5] = aex;
	for (int i = 1, j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = aex;
		val[i + 2] = bex;
		val[i + 4] = cex;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;

	for (int k = 1; k < Nz - 1; k++)//不计算源的部分
	{
		for (int i = 0; i < Nx - 1; i++)
		{
			for (int j = 0; j < Ny - 1; j++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (j == 0 && k != 0)
					rhs[j] = bex*halfgrid_before[k*Nx*Ny + i*Ny + j].ex + cex*halfgrid_before[k*Nx*Ny + i*Ny + j + 1].ex
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bz) / dy - (halfgrid_before[k*Nx*Ny + i*Ny + j].by - halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j].by) / dz);
				else if (k == 0 && j != 0)
					rhs[j] = aex*halfgrid_before[k*Nx*Ny + i*Ny + j - 1].ex + bex*halfgrid_before[k*Nx*Ny + i*Ny + j].ex + cex*halfgrid_before[k*Nx*Ny + i*Ny + j + 1].ex
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bz - halfgrid_before[k*Nx*Ny + i*Ny + j - 1].bz) / dy - (halfgrid_before[k*Nx*Ny + i*Ny + j].by) / dz);
				else if (j == 0 && k == 0)
					rhs[j] = bex*halfgrid_before[k*Nx*Ny + i*Ny + j].ex + cex*halfgrid_before[k*Nx*Ny + i*Ny + j + 1].ex
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bz) / dy - (halfgrid_before[k*Nx*Ny + i*Ny + j].by) / dz);
				else if (j != 0 && k != 0)
					rhs[j] = aex*halfgrid_before[k*Nx*Ny + i*Ny + j - 1].ex + bex*halfgrid_before[k*Nx*Ny + i*Ny + j].ex + cex*halfgrid_before[k*Nx*Ny + i*Ny + j + 1].ex
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bz - halfgrid_before[k*Nx*Ny + i*Ny + j - 1].bz) / dy - (halfgrid_before[k*Nx*Ny + i*Ny + j].by - halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j].by) / dz);
			}//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK) {
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL) {
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK) {
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int j1 = 0; j1 < Ny - 1; j1++)//解的边界修正
			{							
				//----------------------------------------------二阶mur吸收边界条件-----------------------------------------
				if ( i == 0 )//处理上下，还有输入源，三个值应为0的平面
				{
					rhs[j1] = 0.0;
				}
				else if (k == Nz - 2)//右边界，mur吸收边界进行处理
				{
					if (i == Nx - 2)//四条棱边的处理
					{
						rhs[j1] = 0.0;						
					}
					
					else if (j1 == 0)
					{						
						rhs[j1] = halfgrid_before[(k - 1) * Nx*Ny + i * Ny + j1 + 1].ex + ((c*dt - sqrt(dz*dz+dy*dy)) / (c*dt + sqrt(dz*dz + dy*dy)))*(halfgrid_now[(k - 1) * Nx*Ny + i * Ny + j1 + 1].ex - halfgrid_before[k*Nx*Ny + i*Ny + j1].ex);
					}
					else if (j1 == Ny - 2)
					{						
						rhs[j1] =  halfgrid_before[(k - 1) * Nx*Ny + i * Ny + j1 - 1].ex + ((c*dt - sqrt(dz*dz + dy*dy)) / (c*dt + sqrt(dz*dz + dy*dy)))*(halfgrid_now[(k - 1) * Nx*Ny + i * Ny + j1 - 1].ex - halfgrid_before[k*Nx*Ny + i * Ny + j1].ex);
					}
					else
					{	
						rhs[j1] = -halfgrid_beforeX2[(k - 1)*Nx*Ny + i*Ny + j1].ex + (c*dt - dz) / (c*dt + dz)*(halfgrid_now[(k - 1)*Nx*Ny + i*Ny + j1].ex + halfgrid_beforeX2[k*Nx*Ny + i*Ny + j1].ex)
							+ (2 * dz / (c*dt + dz))*(halfgrid_before[k*Nx*Ny + i*Ny + j1].ex + halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j1].ex)
							+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[k*Nx*Ny + (i + 1)*Ny + j1].ex - 2 * halfgrid_before[k*Nx*Ny + i*Ny + j1].ex + halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j1].ex
								+ halfgrid_before[(k - 1)*Nx*Ny + (i+1)*Ny + j1].ex - 2 * halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j1].ex + halfgrid_before[(k - 1)*Nx*Ny + (i - 1)*Ny + j1].ex)
							+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[k*Nx*Ny + i*Ny + j1 + 1].ex - 2 * halfgrid_before[k*Nx*Ny + i*Ny + j1].ex
								+ halfgrid_before[k*Nx*Ny + i*Ny + j1 - 1].ex - halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j1 + 1].ex - 2 * halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j1].ex + halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j1 - 1].ex);				
					}
				}

				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j1].ex = halfgrid_before[k*Nx*Ny + i*Ny + j1].ex;
				halfgrid_before[k*Nx*Ny + i*Ny + j1].ex = halfgrid_now[k*Nx*Ny + i*Ny + j1].ex;
				halfgrid_now[k*Nx*Ny + i*Ny + j1].ex = rhs[j1];
			}
			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
}

void mur2_gsscalc_ey(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	//cout << "源的Ey值---------------->" << halfgrid_now[10 * Ny + 5].ey << endl;
	double aey = (-1 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dz)*(1 / dz));
	double bey = (epsl0 + 2 * (dt / 2)*(dt / 2)*(1 / mur0)*(1 / dz)*(1 / dz));
	double cey = aey;

	int nRet;//GSS函数的返回值
	int N = Nz - 1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nz];
	int ind[3 * (Nz - 1) - 2];
	double val[3 * (Nz - 1) - 2];
	double rhs[Nz - 1];

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
	for (int i = 2, j = 0; i + 2 < 3 * N - 4; i = i + 3)
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
	for (int i = 1, j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = aey;
		val[i + 2] = bey;
		val[i + 4] = cey;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;

	for (int i = 0; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny - 1; j++)
		{
			for (int k = 0; k < Nz - 1; k++)
			{
				if (k == 0 && i != 0)
					rhs[k] = bey*halfgrid_before[k*Nx*Ny + i*Ny + j].ey + cey*halfgrid_before[(k + 1)*Nx*Ny + i*Ny + j].ey
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bx) / dz - (halfgrid_before[k*Nx*Ny + i*Ny + j].bz - halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j].bz) / dx);
				else if (i == 0 && k != 0)
					rhs[k] = aey*halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j].ey + bey*halfgrid_before[k*Nx*Ny + i*Ny + j].ey + cey*halfgrid_before[(k + 1)*Nx*Ny + i*Ny + j].ey
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bx - halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j].bx) / dz - (halfgrid_before[k*Nx*Ny + i*Ny + j].bz) / dx);
				else if (k == 0 && i == 0)
					rhs[k] = bey*halfgrid_before[k*Nx*Ny + i*Ny + j].ey + cey*halfgrid_before[(k + 1)*Nx*Ny + i*Ny + j].ey
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bx) / dz - (halfgrid_before[k*Nx*Ny + i*Ny + j].bz) / dx);
				else if (k != 0 && i != 0)
					rhs[k] = aey*halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j].ey + bey*halfgrid_before[k*Nx*Ny + i*Ny + j].ey + cey*halfgrid_before[(k + 1)*Nx*Ny + i*Ny + j].ey
					+ dt*((halfgrid_before[k*Nx*Ny + i*Ny + j].bx - halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j].bx) / dz - (halfgrid_before[k*Nx*Ny + i*Ny + j].bz - halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j].bz) / dx);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK) {
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL) {
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK) {
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int k1 = 0; k1 < Nz - 1; k1++)
			{
				//修正边界并保存结果				
				if ((i == 0 || i == Nx - 2 ) && (k1 != 0) && (k1 != Nz - 2))
				{
					rhs[k1] = 0.0;
				}	

				else if (k1 == 0)
				{
					rhs[k1] = ((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx)*sin(omega*step*dt);
				}
				 	
				else if (k1 == Nz - 2)//右边界，mur吸收边界进行处理
				{
					if (j == 0 || j == Ny - 2 && i != 0 && i != Nx - 2)
					{
						rhs[k1] = 0.0;
					}									
				    else if (i == 0)
					{
						rhs[k1] = halfgrid_before[(k1 - 1) * Nx*Ny + (i + 1) * Ny + j].ey + ((c*dt - sqrt(dz*dz + dx*dx)) / (c*dt + sqrt(dz*dz + dx*dx)))*(halfgrid_now[(k1 - 1)* Nx*Ny + (i + 1) * Ny + j].ey - halfgrid_before[k1*Nx*Ny + i*Ny + j].ey);
					}
					else if (i == Nx - 2)
					{
						rhs[k1] = halfgrid_before[(k1 - 1) * Nx*Ny + (i - 1) * Ny + j].ey + ((c*dt - sqrt(dz*dz + dx*dx)) / (c*dt + sqrt(dz*dz + dx*dx)))*(halfgrid_now[(k1 - 1) * Nx*Ny + (i - 1) * Ny + j].ey - halfgrid_before[k1*Nx*Ny + i*Ny + j].ey);
					}
					
					else
					{						
						rhs[k1] = -halfgrid_beforeX2[(k1 - 1)*Nx*Ny + i*Ny + j].ex + (c*dt - dz) / (c*dt + dz)*(halfgrid_now[(k1 - 1)*Nx*Ny + i*Ny + j].ex + halfgrid_beforeX2[k1*Nx*Ny + i*Ny + j].ex)
							+ (2 * dz / (c*dt + dz))*(halfgrid_before[k1*Nx*Ny + i*Ny + j].ex + halfgrid_before[(k1 - 1)*Nx*Ny + i*Ny + j].ex)
							+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[k1*Nx*Ny + (i + 1)*Ny + j].ex - 2 * halfgrid_before[k1*Nx*Ny + i*Ny + j].ex + halfgrid_before[k1*Nx*Ny + (i - 1)*Ny + j].ex
								+ halfgrid_before[(k1 - 1)*Nx*Ny + (i + 1)*Ny + j].ex - 2 * halfgrid_before[(k1 - 1)*Nx*Ny + i*Ny + j].ex + halfgrid_before[(k1 - 1)*Nx*Ny + (i - 1)*Ny + j].ex)
							+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[k1*Nx*Ny + i*Ny + j + 1].ex - 2 * halfgrid_before[k1*Nx*Ny + i*Ny + j].ex
								+ halfgrid_before[k1*Nx*Ny + i*Ny + j - 1].ex - halfgrid_before[(k1 - 1)*Nx*Ny + i*Ny + j + 1].ex - 2 * halfgrid_before[(k1 - 1)*Nx*Ny + i*Ny + j].ex + halfgrid_before[(k1- 1)*Nx*Ny + i*Ny + j - 1].ex);
					}
				}
				halfgrid_beforeX2[k1*Nx*Ny + i*Ny + j].ey = halfgrid_before[k1*Nx*Ny + i*Ny + j].ey;
				halfgrid_before[k1*Nx*Ny + i*Ny + j].ey = halfgrid_now[k1*Nx*Ny + i*Ny + j].ey;
				halfgrid_now[k1*Nx*Ny + i*Ny + j].ey = rhs[k1];			
			}

			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
}


void mur2_gsscalc_bz(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	//cout << "源的Bz值---------------->" << halfgrid_now[10 * Ny + 5].bz << endl;
	double ahz = (-1 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dx)*(1 / dx));
	double bhz = (mur0 + 2 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dx)*(1 / dx));
	double chz = ahz;

	int nRet;//GSS函数的返回值
	int N = Nx - 1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nx];
	int ind[3 * (Nx - 1) - 2];
	double val[3 * (Nx - 1) - 2];
	double rhs[Nx - 1];

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
	for (int i = 2, j = 0; i + 2 < 3 * N - 4; i = i + 3)
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
	for (int i = 1, j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = ahz;
		val[i + 2] = bhz;
		val[i + 4] = chz;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;

	for (int k = 1; k < Nz - 1; k++)
	{
		for (int j = 0; j < Ny - 1; j++)
		{
			for (int i = 0; i < Nx - 1; i++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (i == 0)
					rhs[i] = bhz*halfgrid_before[k*Nx*Ny + i*Ny + j].bz + chz*halfgrid_before[k*Nx*Ny + (i + 1)*Ny + j].bz
					+ dt*((halfgrid_now[k*Nx*Ny + i*Ny + j + 1].ex - halfgrid_now[k*Nx*Ny + i*Ny + j].ex) / dy - (halfgrid_now[k*Nx*Ny + (i + 1)*Ny + j].ey - halfgrid_now[k*Nx*Ny + i*Ny + j].ey) / dx);
				else
					rhs[i] = ahz*halfgrid_before[k*Nx*Ny + (i - 1)*Ny + j].bz + bhz*halfgrid_before[k*Nx*Ny + i*Ny + j].bz + chz*halfgrid_before[k*Nx*Ny + (i + 1)*Ny + j].bz
					+ dt*((halfgrid_now[k*Nx*Ny + i*Ny + j + 1].ex - halfgrid_now[k*Nx*Ny + i*Ny + j].ex) / dy - (halfgrid_now[k*Nx*Ny + (i + 1)*Ny + j].ey - halfgrid_now[k*Nx*Ny + i*Ny + j].ey) / dx);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK) {
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL) {
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK) {
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);
			for (int i1 = 0; i1 < Nx - 1; i1++)
			{//保存结果
			 //对结果进行修正

			   if (k == 0)//不改变输入源的值
				{
					rhs[i1] = hm * cos((pi / X)*i1*dx)*cos(omega*step*dt);
				}
							
				halfgrid_beforeX2[k*Nx*Ny + i1*Ny + j].bz = halfgrid_before[k*Nx*Ny + i1*Ny + j].bz;
				halfgrid_before[k*Nx*Ny + i1*Ny + j].bz = halfgrid_now[k*Nx*Ny + i1*Ny + j].bz;
				halfgrid_now[k*Nx*Ny + i1*Ny + j].bz = rhs[i1];

			}

			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}

}

void mur2_gsscalc_bx(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	//cout << "源的Bx值---------------->" << halfgrid_now[10 * Ny + 5].bx << endl;
	double ahx = (-1 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dy)*(1 / dy));
	double bhx = (mur0 + 2 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dy)*(1 / dy));
	double chx = ahx;

	int nRet;//GSS函数的返回值
	int N = Ny - 1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Ny];
	int ind[3 * (Ny - 1) - 2];
	double val[3 * (Ny - 1) - 2];
	double rhs[Ny - 1];

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
	for (int i = 2, j = 0; i + 2 < 3 * N - 4; i = i + 3)
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
	for (int i = 1, j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = ahx;
		val[i + 2] = bhx;
		val[i + 4] = chx;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;

	for (int k = 1; k < Nz - 1; k++)
	{
		for (int i = 0; i < Nx - 1; i++)
		{
			for (int j = 0; j < Ny - 1; j++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (j == 0)
					rhs[j] = bhx*halfgrid_before[k*Nx*Ny + i*Ny + j].bx + chx*halfgrid_before[k*Nx*Ny + i*Ny + j + 1].bx
					+ dt*((halfgrid_now[(k + 1)*Nx*Ny + i*Ny + j].ey - halfgrid_now[k*Nx*Ny + i*Ny + j].ey) / dz - (halfgrid_now[k*Nx*Ny + i*Ny + j + 1].ez - halfgrid_now[k*Nx*Ny + i*Ny + j].ez) / dy);
				else
					rhs[j] = ahx*halfgrid_before[k*Nx*Ny + i*Ny + j - 1].bx + bhx*halfgrid_before[k*Nx*Ny + i*Ny + j].bx + chx*halfgrid_before[k*Nx*Ny + i*Ny + j + 1].bx
					+ dt*((halfgrid_now[(k + 1)*Nx*Ny + i*Ny + j].ey - halfgrid_now[k*Nx*Ny + i*Ny + j].ey) / dz - (halfgrid_now[k*Nx*Ny + i*Ny + j + 1].ez - halfgrid_now[k*Nx*Ny + i*Ny + j].ez) / dy);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK) {
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL) {
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK) {
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int j1 = 0; j1 < Ny - 1; j1++)
			{				
				if (i == 0 || i == Nx - 2 )
				{
					rhs[j1] = 0.0;
				}
				else if (k == 0)
				{
					rhs[j1] = -1 * (X*bate / pi)*hm*sin((pi / X)*i*dx)*sin(omega*step*dt);
				}
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j1].bx = halfgrid_before[k*Nx*Ny + i*Ny + j1].bx;
				halfgrid_before[k*Nx*Ny + i*Ny + j1].bx = halfgrid_now[k*Nx*Ny + i*Ny + j1].bx;
				halfgrid_now[k*Nx*Ny + i*Ny + j1].bx = rhs[j1];				
			}
			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
}

void mur2_gsscalc_by(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now, int step)
{
	//cout << "源的By值---------------->" << halfgrid_now[10 * Ny + 5].by<< endl;
	//inject_field(halfgrid_before, halfgrid_now, step);
	double ahy = (-1 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dz)*(1 / dz));
	double bhy = (mur0 + 2 * (dt / 2)*(dt / 2)*(1 / epsl0)*(1 / dz)*(1 / dz));
	double chy = ahy;

	int nRet;//GSS函数的返回值
	int N = Nz - 1;
	int nnz = 3 * N - 2;
	int nRow = N;
	int nCol = N;

	int ptr[Nz];
	int ind[3 * (Nz - 1) - 2];
	double val[3 * (Nz - 1) - 2];
	double rhs[Nz - 1];

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
	for (int i = 2, j = 0; i + 2 < 3 * N - 4; i = i + 3)
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
	for (int i = 1, j = 2; i + 3 < 3 * N - 2; i = i + 3)
	{
		val[i] = ahy;
		val[i + 2] = bhy;
		val[i + 4] = chy;
		j++;
	}
	//-----------------------------------------------------------rhs数组初始化
	for (int i = 0; i < N; i++)
		rhs[i] = 0.0;

	for (int i = 0; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny - 1; j++)
		{
			for (int k = 0; k < Nz - 1; k++)
			{
				//首先处理矩阵方程的右端项rhs数组
				if (k == 0)
					rhs[k] = bhy*halfgrid_before[k*Nx*Ny + i*Ny + j].by + chy*halfgrid_before[(k + 1)*Nx*Ny + i*Ny + j].by
					+ dt*((halfgrid_now[k*Nx*Ny + (i + 1)*Ny + j].ez - halfgrid_now[k*Nx*Ny + i*Ny + j].ez) / dx - (halfgrid_now[(k + 1)*Nx*Ny + i*Ny + j].ex - halfgrid_now[k*Nx*Ny + i*Ny + j].ex) / dz);
				else
					rhs[k] = ahy*halfgrid_before[(k - 1)*Nx*Ny + i*Ny + j].by + bhy*halfgrid_before[k*Nx*Ny + i*Ny + j].by + chy*halfgrid_before[(k + 1)*Nx*Ny + i*Ny + j].by
					+ dt*((halfgrid_now[k*Nx*Ny + (i + 1)*Ny + j].ez - halfgrid_now[k*Nx*Ny + i*Ny + j].ez) / dx - (halfgrid_now[(k + 1)*Nx*Ny + i*Ny + j].ex - halfgrid_now[k*Nx*Ny + i*Ny + j].ex) / dz);
			}
			//GSS求解
			nRet = GSS_init_ld(nRow, nCol, ptr, ind, val, type, setting);
			if (nRet != GRUS_OK) {
				printf("\tERROR at init GSS solver. ERROR CODE:%d\r\n", nRet);
				return;
			}

			hSolver = GSS_symbol_ld(nRow, nCol, ptr, ind, val);
			if (hSolver == NULL) {
				printf("\tERROR at SYMBOLIC ANALYSIS.\r\n");
				exit(0);
			}

			nRet = GSS_numeric_ld(nRow, nCol, ptr, ind, val, hSolver);
			if (nRet != GRUS_OK) {
				printf("\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n", nRet);
				hSolver = NULL;		//必须设置为NULL,GSS已自动释放内存
				exit(0);
			}

			GSS_solve_ld(hSolver, nRow, nCol, ptr, ind, val, rhs);

			for (int k1 = 0; k1 < Nz - 1; k1++)
			{//保存结果				
				if (k1 == 0 || j == 0 || j == Ny - 2)
				{
					rhs[k1] = 0.0;
				}

				halfgrid_beforeX2[k1*Nx*Ny + i*Ny + j].by = halfgrid_before[k1*Nx*Ny + i*Ny + j].by;
				halfgrid_before[k1*Nx*Ny + i*Ny + j].by = halfgrid_now[k1*Nx*Ny + i*Ny + j].by;
				halfgrid_now[k1*Nx*Ny + i*Ny + j].by = rhs[k1];				
			}
			if (hSolver != NULL)
				GSS_clear_ld(hSolver);
		}
	}
}

#endif