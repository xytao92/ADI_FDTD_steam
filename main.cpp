﻿//writen by liuzhichao 
//2016/10/8
#include "definer.h"
#include "grid.h"
#include "functions.h"
#include"GSS_ADI_matel.h"
#include"GSS_ADI_mur.h"

int main()
{
	clock_t start_time = clock();
	initGrid(halfgrid_beforeX2,halfgrid_before, halfgrid_now);//对网格数据进行初始化
	int choice;
	
	while (TRUE)
	{
		cout << "请输入吸收边界类型:\n" << endl;
		cout << "金属吸收边界请按--1\n" << endl;
		cout << "Mur吸收边界请按--2(维护中)\n" << endl;
		cin >> choice;
		if (choice == 1)
		{
			adi_fdtd_leapforg_matel_GSS(halfgrid_before, halfgrid_now);//进行计算
		break;
	    }
		else if (choice == 2)
		{
			adi_fdtd_leapforg_mur_GSS(halfgrid_beforeX2, halfgrid_before, halfgrid_now);//进行计算
			break;
		}
		else
		{
			cout << "抱歉，无此选项，请重新输入。" << endl;
		}
	}
	
	free(halfgrid_beforeX2, halfgrid_before, halfgrid_now);//保存结果,并释放空间

	clock_t end_time = clock();
	cout << "Running time: " << static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	system("pause");
	return 0;
}