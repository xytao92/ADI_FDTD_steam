//writen by liuzhichao 
//2016/10/8

#ifndef _FUNCTIONS_
#define _FUNCTIONS_
#include "definer.h"
#include "grid.h"

//============================初始化网格==========================//
void initGrid(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now)//网格的初始化
{
	for (int k = 0; k<Nz; k++)
		for (int i = 0; i<Nx; i++)
			for (int j = 0; j<Ny; j++)
			{
				halfgrid_now[k*Nx*Ny + i*Ny + j].ex = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].ey = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].ez = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].bx = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].by = 0.0;
				halfgrid_now[k*Nx*Ny + i*Ny + j].bz = 0.0;

				halfgrid_before[k*Nx*Ny + i*Ny + j].ex = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].ey = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].ez = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].bx = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].by = 0.0;
				halfgrid_before[k*Nx*Ny + i*Ny + j].bz = 0.0;

				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].ex = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].ey = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].ez = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].bx = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].by = 0.0;
				halfgrid_beforeX2[k*Nx*Ny + i*Ny + j].bz = 0;

				/*grid_result[k*Nx*Ny + i*Ny + j].ex = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].ey = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].ez = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].bx = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].by = 0.0;
				grid_result[k*Nx*Ny + i*Ny + j].bz = 0.0;*/
			}

}

//==========================初始化输入源==========================//
void init_source(Grid* halfgrid_now)//计算激励源//
{
	int k =0;//加点源的位置//
	int step = 0;
	for (int i = 1; i < Nx - 2; i++)
	{
		for (int j = 1; j < Ny - 2; j++)
		{
			if (step*dt < 2 * T)
			{
				//加一个完整的TE10,模式

				//电场
				double rising_edge = (step*dt) / (2 * T);

				/*halfgrid_now[i*Ny + j].ez = rising_edge * 0.0;*/

				halfgrid_now[i*Ny + j].ex = rising_edge * 0.0;

				double temp_ey = sin((pi / X)*i*dx);
				halfgrid_now[i*Ny + j].ey = rising_edge * temp_ey * sin(omega*step*dt);

				//磁场
				/*double temp_bz = hm*cos((pi / X)*i*dx);
				halfgrid_now[i*Ny + j].bz = rising_edge * temp_bz * cos(omega*step*dt);*/

				/*double temp_bx = (X*bate / pi)*hm*sin((pi / X)*i*dx);*/

				/*halfgrid_now[i*Ny + j].bx = -1 * rising_edge * temp_bx * sin(omega*step*dt*1.5);*/

				/*halfgrid_now[i*Ny + j].bx = rising_edge * 0.0;*/

				/*halfgrid_now[i*Ny + j].by = rising_edge * 0.0;*/

			}
			else 
			{
				//电场
				/*halfgrid_now[i*Ny + j].ez = 0.0;*/

				halfgrid_now[i*Ny + j].ex = 0.0;

				double temp_ey =  sin((pi / X)*i*dx);
				halfgrid_now[i*Ny + j].ey = temp_ey * sin(omega*step*dt);

				/*double temp_ey = ((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx);
				halfgrid_now[i*Ny + j].ey = temp_ey * sin(omega*step*dt); *///-1 * pi*sin(pi*i*dx / X)*sin(2 * pi*freq*step*dt) / X;//((omega*mur0*X) / pi) * hm * sin((pi / X)*i*dx) * sin(omega*step*dt);
				
		       //磁场
				/*double temp_bz = hm*cos((pi / X)*i*dx);
				halfgrid_now[i*Ny + j].bz = temp_bz * cos(omega*step*dt);*/

				/*double temp_bx = -1 * (X*bate / pi)*hm*sin((pi / X)*i*dx); */
				/*halfgrid_now[i*Ny + j].bx = temp_bx * sin(omega*step*dt*1.5);*/

				/*halfgrid_now[i*Ny + j].bx =0.0;*/

				/*halfgrid_now[i*Ny + j].by = 0.0;*/
			
			}
			
		}
	}
}

//=============================保存结果===========================//
void save_result(ofstream file, Grid* halfgrid_now, int step)
{
	int result_z = 2;
	int result_x = 10;
	int result_y = 5;

	file << step << '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ex
		<< '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ey
		<< '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].ez << '\t';

	file << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].bx
		<< '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].by
		<< '\t' << halfgrid_now[result_z * Nx*Ny + result_x * Ny + result_y].bz << '\t';

	file << '\n';

}

//=========================输出横截面数据=========================//
void get_plat(ofstream file2,Grid* halfgrid_now)//获取一个截面的数据
{
	for (int i = 0; i < Nx - 1; i++)
	{
		for (int k = 0; k < Nz - 1; k++)
		{
			file2 << k << '\t' << i << '\t'
				<< halfgrid_now[k*Nx*Ny + i*Ny + 5].ey << '\t'
				<< halfgrid_now[k*Nx*Ny + i*Ny + 5].bx << '\t'
				<< halfgrid_now[k*Nx*Ny + i*Ny + 5].bz << endl;
		}
	}

}

//============================释放空间===========================//
void free(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now)
{
	delete halfgrid_before;
	delete halfgrid_now;
	delete halfgrid_beforeX2;
	//delete grid_result;
}

//=======================计算CFL稳定性条件=======================//
double CFL_calc()
{
	double result = 0.0;
	double v = 1 / (sqrt(epsl0*mur0));
	double g_temp = sqrt((1 / dx)*(1 / dx) + (1 / dy)*(1 / dy) + (1 / dz)*(1 / dz));
	result = 1 / (v*g_temp);
	return result;
}
double fc_calc()
{
	double fc = 1 / (2 * X * sqrt(epsl0*mur0));
	return fc;
}

//===========================输出理论值=========================//
int theor_val_gen()
{
	clock_t start_time = clock();
	system("mkdir result\\theor_val");
	double ex = 0.0;
	double ey = 0.0;
	double ez = 0.0;
	double hx = 0.0;
	double hy = 0.0;
	double hz = 0.0;
	double temp0 = 0.0;
	double temp1 = 0.0;
	double temp2 = 0.0;

	int step = 0;

	int result_z = 10;
	int result_x = 10;
	int result_y = 5;//1270

	ofstream file(theor_val_filepath);//用于保存结果

	while (step < STEPS)
	{
		if (step*dt < 2 * T)
		{
			//=============================电场=============================//
			ex = ((step*dt) / (2 * T))*0.0;//ex

			//double temp0 = ((omega*mur0*X) / pi) * hm * sin((pi / X)*result_x*dx);
			double temp0 =  hm * sin((pi / X)*result_x*dx);
			double temp_ey_sin = omega*step*dt - bate*result_z;
			double temp_ey_font = (step*dt) / (2 * T);

			ey = temp_ey_font * temp0 * sin(temp_ey_sin);//ey


			ez = ((step*dt) / (2 * T))*0.0;//ez

			//=============================磁场=============================//

			//double temp1 = (X*bate / pi) * sin((pi / X) * result_x * dx);
			double temp1 = (bate / (omega*mur0))*sin((pi / X) * result_x * dx);
			double temp_hx_sin = omega*step*dt - bate*result_z;
			double temp_hx_font = -1 * (step*dt) / (2 * T) * hm;

			hx = temp_hx_font * temp1 * sin(temp_hx_sin);//hx

			hy = ((step*dt) / (2 * T))*0.0;//hy

			double temp4 = (pi / X) * result_x*dx;
			double temp2 = (pi / (omega*mur0*X))*hm * cos(temp4);
			double temp3 = omega*step*dt - bate*result_z;
			double temp_hz_font = (step*dt) / (2 * T);

			hz = temp_hz_font * temp2 * cos( temp3 );//hz
		}
		else
		{
			//=============================电场=============================//
			ex = 0.0;//ex
			double temp_ey_sin_1 = (pi / X)*result_x*dx;
			//double temp0 = ((omega*mur0*X) / pi) * hm * sin(temp_ey_sin_1);
			double temp0 = hm * sin(temp_ey_sin_1);
			double temp_ey_sin_2 = omega*step*dt - bate*result_z;
			ey = temp0*sin(temp_ey_sin_2);//ey
			ez = 0.0;//ez
			//=============================磁场=============================//
			double temp_hx_1 = (pi / X) * result_x * dx;
			//double temp1 = (X*bate / pi)* hm* sin(temp_hx_1);
			double temp1 = (bate / (omega*mur0))* hm* sin(temp_hx_1);
			double temp_hx_2 = omega*step*dt - bate*result_z;
			hx = -1 * temp1 * sin(temp_hx_2);//hx
			//由于bate会在不同的频率下产生误差会造成相当于输入了两个源
			hy = 0.0;//hy
			double temp2 = (pi / (omega*mur0*X))* hm * cos((pi / X) * result_x*dx);
			double temp_hz_cos = omega*step*dt - bate*result_z;
			hz = temp2 * cos(temp_hz_cos);//hz
		}
		file << step << '\t' << ex << '\t' << ey << '\t' << ez << '\t' << hx << '\t' << hy << '\t' << hz << endl;
		cout << "STEP---" << step << "  has finished" << endl;
		step++;
	}
	
	clock_t end_time = clock();
	cout << "理论值计算所用时间 : " << static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC << " s\n" << endl;
	return 0;
}











//以下函数没有修改坐标系，按照老的坐标系来------------------------！！！！！！！

//void adi_fdtd_leapforg_matel(Grid* halfgrid_before,Grid* halfgrid_now)
//
//{
//	/*system("mkdir result");*/
//	ofstream file("result\\leapforg_ADI_FDTD_steam0.1_matel03.txt");//用于保存结果
//	int i1 = 0;
//	int j1 = 0;
//	int k1 = 0;
//	//*******计算TE10模******//注意边界条件的问题，不处理周围的四个面
//	//PART1---- 计算电场//
//
//	int step = 0;//计算时间步长
//	while(step< STEPS)
//	{
//		if (step == 50)
//		{
//			cout << "It's Time, step =50!" << endl;
//		}
//
//		inject_field(halfgrid_before,halfgrid_now,step);
//		for(int i = 0; i<Nx-1;i++)
//			for(int j =1; j<Ny-2; j++)
//				for (int k = 1; k<Nz-2; k++)
//				{
//			i1 = i;
//			j1 = j;
//			k1 = k;
//			if (i == 100 && j == 100 && k == 50)
//			{
//				cout << "Now" << endl;
//			}
//			if (i == 0)//左边界，信号源处
//			{
//				i1 = 1;
//			}
//			
//			//ex
//			double ty0 = ( halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[i*Ny*Nz + (j1-1)*Nz + k].bz )/dy;
//			double tz0 = ( halfgrid_before[i*Ny*Nz + j*Nz + k].by-halfgrid_before[i*Ny*Nz + j*Nz + k1-1].by )/dz;
//			halfgrid_now[ i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i*Ny*Nz + j*Nz + k].ex + dt*(1/epsl_x)*( ty0 - tz0 );
//			//ey
//			double tx1 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + j*Nz + k1-1].bx )/dz;
//			double tz1 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].bz -halfgrid_before[(i1-1)*Ny*Nz + j*Nz + k].bz )/dx;
//			halfgrid_now[ i*Ny*Nz + j*Nz + k].ey = halfgrid_before[ i*Ny*Nz + j*Nz + k].ey + dt*(1/epsl_y)*( tx1 - tz1 );
//			//ez
//			double ty2 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[(i1-1)*Ny*Nz + j*Nz + k ].by )/dx;
//			double tx2 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k].bx-halfgrid_before[i*Ny*Nz + (j1-1)*Nz + k].bx )/dy;
//			halfgrid_now[ i*Ny*Nz + j*Nz + k].ez = halfgrid_before[ i*Ny*Nz + j*Nz + k].ez + dt*(1/epsl_z)*( ty2 - tx2 );
//
//			halfgrid_before[ i*Ny*Nz + j*Nz + k].ex = halfgrid_now[ i*Ny*Nz + j*Nz + k].ex;//保存为前一步结果
//			halfgrid_before[ i*Ny*Nz + j*Nz + k].ey =  halfgrid_now[i*Ny*Nz + j*Nz + k].ey;
//			halfgrid_before[i*Ny*Nz + j*Nz + k].ez  = halfgrid_now[i*Ny*Nz + j*Nz + k].ez;
//
//			//bx
//			double ty3 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dz;
//			double tz3 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dy;
//			halfgrid_now[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[i*Ny*Nz + j*Nz + k].bx + dt*(1 / mur_x)*(ty3 - tz3);
//			//by
//			double tz4 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dx;
//			double tx4 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dz;
//			halfgrid_now[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by + dt*(1 / mur_y)*(tz4 - tx4);
//			//bz
//			double tx5 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dy;
//			double ty5 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dx;
//			halfgrid_now[i*Ny*Nz + j*Nz + k].bz = halfgrid_before[i*Ny*Nz + j*Nz + k].bz + dt*(1 / mur_z)*(tx5 - ty5);
//
//			halfgrid_before[i*Ny*Nz + j*Nz + k].bx = halfgrid_now[i*Ny*Nz + j*Nz + k].bx;//保存为前一步结果
//			halfgrid_before[i*Ny*Nz + j*Nz + k].by = halfgrid_now[i*Ny*Nz + j*Nz + k].by;
//			halfgrid_before[i*Ny*Nz + j*Nz + k].bz = halfgrid_now[i*Ny*Nz + j*Nz + k].bz;
//			
			
//PART2----计算磁场//
	//for(int i = 0; i<Nx-1;i++)
	//	for(int j =0; j<Ny-1; j++)
	//		for (int k = 0; k < Nz-1; k++)
	//		{
	//			//bx
	//			double ty3 =  ( halfgrid_before[i*Ny*Nz + j*Nz + k+1].ey - halfgrid_before[ i*Ny*Nz + j*Nz + k].ey )/dz;
	//			double tz3 =  ( halfgrid_before[ i*Ny*Nz + (j+1)*Nz + k].ez-halfgrid_before[ i*Ny*Nz +  j*Nz + k ].ez  )/dy;
	//			halfgrid_now[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[ i*Ny*Nz + j*Nz + k].bx + dt*(1/ mur_x )*( ty3 - tz3 );
	//			//by
	//			double tz4 =  ( halfgrid_before[ (i+1)*Ny*Nz + j*Nz + k].ez - halfgrid_before[ i*Ny*Nz +  j*Nz + k].ez )/dx;
	//			double tx4 =  (halfgrid_before[i*Ny*Nz + j*Nz + k+1].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex )/dz;
	//			halfgrid_now[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by + dt*(1/ mur_y )*( tz4 - tx4 );
	//			//bz
	//			double tx5 =  ( halfgrid_before[ i*Ny*Nz + (j+1)*Nz + k].ex - halfgrid_before[ i*Ny*Nz + j*Nz + k ].ex )/dy;
	//			double ty5 =  (halfgrid_before[ (i+1)*Ny*Nz + j*Nz + k].ey - halfgrid_before[ i*Ny*Nz + j*Nz + k].ey )/dx;
	//			halfgrid_now[ i*Ny*Nz + j*Nz + k].bz = halfgrid_before[ i*Ny*Nz + j*Nz + k].bz + dt*(1/ mur_z )*( tx5 - ty5 );			

	//			halfgrid_before[i*Ny*Nz + j*Nz + k].bx = halfgrid_now[i*Ny*Nz + j*Nz + k].bx;//保存为前一步结果
 //               halfgrid_before[ i*Ny*Nz + j*Nz + k].by =  halfgrid_now[i*Ny*Nz + j*Nz + k].by;
	//			halfgrid_before[ i*Ny*Nz + j*Nz + k].bz  = halfgrid_now[ i*Ny*Nz + j*Nz + k].bz;
			//}
								
	////若取平均，则网格点处的电磁场计算如下
	//
	//for(int i = 0; i<Nx-1;i++)
	//	for(int j =0; j<Ny-1; j++)
	//		for (int k = 0; k < Nz-1; k++)
	//		{
	//			  if(i == 0 )//左边界，信号源处
	//			 {
	//				
	//			 }
	//		           
	//		     if(i== Nx - 1 )//右边界，吸收边界
	//			 {
	//				
	//			 }
	//			 if(k == Nz - 1|| k==0 || j==Ny - 1 || j==0 )//四个平面
	//			 {
	//				 
	//			 }
	//			 else
	//			 {
	//			 
	//			    grid_result[i*Ny*Nz + j*Nz + k].ex=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k].ex + halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+2].ex)/4;
	//				grid_result[i*Ny*Nz + j*Nz + k].ey=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k].ey +  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k].ey + halfgrid_now[ i*Ny*Nz +( j+1) *Nz + k+2].ey+  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k+2].ey )/4;
	//				grid_result[i*Ny*Nz + j*Nz + k].ez= ( halfgrid_now[ i*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ i*Ny*Nz + (j+2)*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + (j+2)*Nz + k+1].ez)/4;
	//				grid_result[i*Ny*Nz + j*Nz + k].bx=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k+1].bx + halfgrid_now[ (i+2)*Ny*Nz + (j+1)*Nz + k+2].bx  )/2;
	//				grid_result[i*Ny*Nz + j*Nz + k].by=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].by +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+1].by  )/2;
	//				grid_result[i*Ny*Nz + j*Nz + k].bz=(  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k ].bz +  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k+2].bz)/2;
	//			 }//else
	//		}	
//	int result_x = 100;
//	int result_y = 100;
//	int result_z = 50;
//	file << step << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ex << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ey << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ez << '\t';
//	file << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bx << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].by << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bz << '\t';
//	file << '\n';
//	//save_result(halfgrid_now, step);
//	step++;
//	cout<<"Step--- "<<step<<" ---has finished."<<endl;
// }//while
//
//}//函数结尾
//
//void adi_fdtd_leapforg_mur(Grid* halfgrid_beforeX2,Grid* halfgrid_before, Grid* halfgrid_now)
//{
//	/*system("mkdir result");*/
//	ofstream file("result\\leapforg_ADI_FDTD_steam0.1_mur02.txt");//用于保存结果
//	//*******计算TE10模******//注意边界条件的问题,直接全部使用二阶mur吸收边界还是会出现问题，因为有棱边的存在，现在采用在棱边处使用一阶吸收边界条件，避免使用棱边
//
//	int step = 0;//计算时间步长
//	while (step< STEPS)
//	{
//		inject_field(halfgrid_before,halfgrid_now, step);
//		//PART1---- 计算电场//
//		for (int i = 0; i<Nx - 1; i++)
//			for (int j = 1; j<Ny - 2; j++)
//				for (int k = 1; k<Nz - 2; k++)
//				{
//			if (i == 0)//左边界，信号源处,一阶mur吸收边界处理棱边			
//			{
//				/*if (j == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i+1) * Ny*Nz + (j+1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i+1)* Ny*Nz + (j+1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//				}
//				else if (j == Ny - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i+1) * Ny*Nz + (j-1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i+1) * Ny*Nz + (j-1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//				}
//				else if (k == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k+1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i + 1) * Ny*Nz + j * Nz + k+1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//				}
//				else if (k == Nz - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i + 1) * Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//				}*/
//			  //else
//				//{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[(i + 1)*Ny*Nz + j*Nz + k].ey + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
//					+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey
//					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k].ey)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
//					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey - halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ey);
//					
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[(i + 1)*Ny*Nz + j*Nz + k].ez + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i + 1)*Ny*Nz + j*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
//					+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez
//					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k - 1].ez)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
//					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ez);
//				//}
//
//			}
//			else if (i == Nx - 1)//右边界，信号源处
//			{
//					/*if (j == 0)
//					{
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j + 1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1)* Ny*Nz + (j + 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//					}
//					else if (j == Ny - 1)
//					{
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//					}
//					else if (k == 0)
//					{
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k + 1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k + 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//					}
//					else if (k == Nz - 1)
//					{
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dx) / (c*dt + sqrt(2)*dx))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//					}
//					else
//					{*/
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//						halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[(i - 1)*Ny*Nz + j*Nz + k].ey + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i - 1)*Ny*Nz + j*Nz + k].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
//						+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey)
//						+ (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey
//					    + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ey)
//					    + (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
//					    + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k - 1].ey);
//
//				        halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[(i - 1)*Ny*Nz + j*Nz + k].ez + (c*dt - dx) / (c*dt + dx)*(halfgrid_now[(i - 1)*Ny*Nz + j*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
//				    	+ (2 * dx / (c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez)
//					    + (dx*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez
//					    + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ez)
//					    + (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
//					    + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k - 1].ez);
//
//					//}
//															
//			}
//			/*else if (j == 0)
//			{
//				if (i == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i + 1) * Ny*Nz + (j + 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i + 1)* Ny*Nz + (j + 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//				}
//				else if (i == Nx - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j + 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j + 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//				}
//				else if (k == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j + 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j + 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else if (k == Nz - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j+1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[i * Ny*Nz + (j+1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else
//				{				
//			    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//
//				halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + (j + 1)*Nz + k].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
//					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex)
//					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + (j + 1)*Nz + k].ex)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ex);
//
//				halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[i*Ny*Nz + (j + 1)*Nz + k].ez + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j + 1)*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
//					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez)
//					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez
//					+ halfgrid_before[(i + 1)*Ny*Nz + (j + 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + (j + 1)*Nz + k].ez)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
//					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ez);
//
//				}
//				
//			}
//			else if (j == Ny - 1)
//			{
//				if (i == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i + 1) * Ny*Nz + (j - 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i + 1)* Ny*Nz + (j - 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//				}
//				else if (i == Nx - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j - 1) * Nz + k].ez - halfgrid_before[i * Ny*Nz + j * Nz + k].ez);
//				}
//				else if (k == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[(i - 1) * Ny*Nz + (j - 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else if (k == Nz - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dy) / (c*dt + sqrt(2)*dy))*(halfgrid_now[i * Ny*Nz + (j - 1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else
//				{
//				    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//
//				    halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + (j - 1)*Nz + k].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j - 1)*Nz + k].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
//					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex)
//					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ex)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ex);
//
//				    halfgrid_now[i*Ny*Nz + j*Nz + k].ez = -halfgrid_beforeX2[i*Ny*Nz + (j - 1)*Nz + k].ez + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + (j - 1)*Nz + k].ez + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez)
//					+ (2 * dy / (c*dt + dy))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez)
//					+ (dy*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dy))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ez
//					+ halfgrid_before[(i + 1)*Ny*Nz + (j - 1)*Nz + k].ez - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez + halfgrid_before[(i - 1)*Ny*Nz + (j - 1)*Nz + k].ez)
//					+ (dx*(c*dt)*(c*dt)) / (2 * dz*dz*(c*dt + dx))*(halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ez
//					+ halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ez - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ez - 2 * halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ez + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ez);
//
//				}
//			}
//
//			else if (k == 0)
//			{
//				if (i == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k + 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i + 1)* Ny*Nz + j * Nz + k + 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//				}
//				else if (i == Nx - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k + 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k+1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//				}
//				else if (j == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j + 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j + 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else if (j == Ny - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k + 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j - 1) * Nz + k + 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//
//				    halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k + 1].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
//					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ex);
//
//				    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k + 1].ey + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k + 1].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
//					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey
//					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey + halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
//					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k + 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k + 1].ey);
//				}
//
//				
//
//			}
//			else if (k == Nz - 1)
//			{
//				if (i == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i + 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i + 1)* Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//				}
//				else if (i == Nx - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[(i - 1) * Ny*Nz + j * Nz + k - 1].ey + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[(i - 1) * Ny*Nz + j * Nz + k - 1].ey - halfgrid_before[i * Ny*Nz + j * Nz + k].ey);
//				}
//				else if (j == 0)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j + 1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j + 1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else if (j == Ny - 1)
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ey = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i * Ny*Nz + (j - 1) * Nz + k - 1].ex + ((c*dt - sqrt(2)*dz) / (c*dt + sqrt(2)*dz))*(halfgrid_now[i * Ny*Nz + (j - 1) * Nz + k - 1].ex - halfgrid_before[i * Ny*Nz + j * Nz + k].ex);
//				}
//				else
//				{
//					halfgrid_now[i*Ny*Nz + j*Nz + k].ez = 0;
//
//				    halfgrid_now[i*Ny*Nz + j*Nz + k].ex = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k - 1].ex + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k - 1].ex + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex)
//					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ex
//					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ex - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ex + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ex);
//
//				    halfgrid_now[i*Ny*Nz + j*Nz + k].ey = -halfgrid_beforeX2[i*Ny*Nz + j*Nz + k - 1].ey + (c*dt - dy) / (c*dt + dy)*(halfgrid_now[i*Ny*Nz + j*Nz + k - 1].ey + halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey)
//					+ (2 * dz / (c*dt + dz))*(halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dx*dx*(c*dt + dz))*(halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey + halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].ey
//					+ halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k - 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey + halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey)
//					+ (dz*(c*dt)*(c*dt)) / (2 * dy*dy*(c*dt + dz))*(halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k].ey
//					+ halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].ey - halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k - 1].ey - 2 * halfgrid_before[i*Ny*Nz + j*Nz + k - 1].ey + halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k - 1].ey);
//
//				}
//				
//			}
//			*/
//
//			else
//			{
//
//				//ex
//				double ty0 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bz) / dy;
//				double tz0 = (halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].by) / dz;
//				halfgrid_now[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i*Ny*Nz + j*Nz + k].ex + dt*(1 / epsl_x)*(ty0 - tz0);
//				//ey
//				double tx1 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + j*Nz + k - 1].bx) / dz;
//				double tz1 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bz - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].bz) / dx;
//				halfgrid_now[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[i*Ny*Nz + j*Nz + k].ey + dt*(1 / epsl_y)*(tx1 - tz1);
//				//ez
//				double ty2 = (halfgrid_before[i*Ny*Nz + j*Nz + k].by - halfgrid_before[(i - 1)*Ny*Nz + j*Nz + k].by) / dx;
//				double tx2 = (halfgrid_before[i*Ny*Nz + j*Nz + k].bx - halfgrid_before[i*Ny*Nz + (j - 1)*Nz + k].bx) / dy;
//				halfgrid_now[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[i*Ny*Nz + j*Nz + k].ez + dt*(1 / epsl_z)*(ty2 - tx2);
//
//			}//else
//			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ex = halfgrid_before[i*Ny*Nz + j*Nz + k].ex;//保存上上步的结果
//			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ey = halfgrid_before[i*Ny*Nz + j*Nz + k].ey;
//			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].ez = halfgrid_before[i*Ny*Nz + j*Nz + k].ez;
//
//
//			halfgrid_before[i*Ny*Nz + j*Nz + k].ex = halfgrid_now[i*Ny*Nz + j*Nz + k].ex;//保存前一步结果
//			halfgrid_before[i*Ny*Nz + j*Nz + k].ey = halfgrid_now[i*Ny*Nz + j*Nz + k].ey;
//			halfgrid_before[i*Ny*Nz + j*Nz + k].ez = halfgrid_now[i*Ny*Nz + j*Nz + k].ez;
//				}
//
//		//PART2----计算磁场//
//		for (int i = 0; i<Nx - 1; i++)
//			for (int j = 0; j<Ny - 1; j++)
//				for (int k = 0; k < Nz - 1; k++)
//				{
//
//			//bx
//			double ty3 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dz;
//			double tz3 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dy;
//			halfgrid_now[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[i*Ny*Nz + j*Nz + k].bx + dt*(1 / mur_x)*(ty3 - tz3);
//			//by
//			double tz4 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ez - halfgrid_before[i*Ny*Nz + j*Nz + k].ez) / dx;
//			double tx4 = (halfgrid_before[i*Ny*Nz + j*Nz + k + 1].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dz;
//			halfgrid_now[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by + dt*(1 / mur_y)*(tz4 - tx4);
//			//bz
//			double tx5 = (halfgrid_before[i*Ny*Nz + (j + 1)*Nz + k].ex - halfgrid_before[i*Ny*Nz + j*Nz + k].ex) / dy;
//			double ty5 = (halfgrid_before[(i + 1)*Ny*Nz + j*Nz + k].ey - halfgrid_before[i*Ny*Nz + j*Nz + k].ey) / dx;
//			halfgrid_now[i*Ny*Nz + j*Nz + k].bz = halfgrid_before[i*Ny*Nz + j*Nz + k].bz + dt*(1 / mur_z)*(tx5 - ty5);
//
//			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].bx = halfgrid_before[i*Ny*Nz + j*Nz + k].bx;//保存上上步的结果
//			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].by = halfgrid_before[i*Ny*Nz + j*Nz + k].by;
//			halfgrid_beforeX2[i*Ny*Nz + j*Nz + k].bz = halfgrid_before[i*Ny*Nz + j*Nz + k].bz;
//
//			halfgrid_before[i*Ny*Nz + j*Nz + k].bx = halfgrid_now[i*Ny*Nz + j*Nz + k].bx;//保存为前一步结果
//			halfgrid_before[i*Ny*Nz + j*Nz + k].by = halfgrid_now[i*Ny*Nz + j*Nz + k].by;
//			halfgrid_before[i*Ny*Nz + j*Nz + k].bz = halfgrid_now[i*Ny*Nz + j*Nz + k].bz;
//				}
//
//		////若取平均，则网格点处的电磁场计算如下
//		//
//		//for(int i = 0; i<Nx;i++)
//		//	for(int j =0; j<Ny; j++)
//		//		for (int k = 0; k < Nz; k++)
//		//		{
//		//			  if(i == 0 )//左边界，信号源处
//		//			 {
//		//				
//		//			 }
//		//		           
//		//		     if(i== Nx - 1 )//右边界，吸收边界
//		//			 {
//		//				
//		//			 }
//		//			 if(k == Nz - 1|| k==0 || j==Ny - 1 || j==0 )//四个平面
//		//			 {
//		//				 
//		//			 }
//		//			 else
//		//			 {
//		//			 
//		//			    grid_result[i*Ny*Nz + j*Nz + k].ex=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k].ex + halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k].ex +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+2].ex)/4;
//		//				grid_result[i*Ny*Nz + j*Nz + k].ey=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k].ey +  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k].ey + halfgrid_now[ i*Ny*Nz +( j+1) *Nz + k+2].ey+  halfgrid_now[ (i+2)*Ny*Nz +( j+1) *Nz + k+2].ey )/4;
//		//				grid_result[i*Ny*Nz + j*Nz + k].ez= ( halfgrid_now[ i*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ i*Ny*Nz + (j+2)*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + j*Nz + k+1].ez + halfgrid_now[ (i+2)*Ny*Nz + (j+2)*Nz + k+1].ez)/4;
//		//				grid_result[i*Ny*Nz + j*Nz + k].bx=(  halfgrid_now[ i*Ny*Nz + (j+1)*Nz + k+1].bx + halfgrid_now[ (i+2)*Ny*Nz + (j+1)*Nz + k+2].bx  )/2;
//		//				grid_result[i*Ny*Nz + j*Nz + k].by=(  halfgrid_now[ (i+1)*Ny*Nz + j*Nz + k+1].by +  halfgrid_now[ (i+1)*Ny*Nz + (j+2)*Nz + k+1].by  )/2;
//		//				grid_result[i*Ny*Nz + j*Nz + k].bz=(  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k ].bz +  halfgrid_now[ (i+1)*Ny*Nz + (j+1)*Nz + k+2].bz)/2;
//		//			 }//else
//		//		}	
//		int result_x = 100; 
//		int result_y = 100;
//		int result_z = 50;
//		file << step << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ex << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ey << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].ez << '\t';
//		file << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bx << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].by << '\t' << halfgrid_now[result_x * Ny*Nz + result_y * Nz + result_z].bz << '\t';
//		file << '\n';
//		//save_result(halfgrid_now, step);
//		step++;
//		cout << "Step--- " << step << " ---has finished." << endl;
//	}//while
//
//}//函数结尾
//

#endif