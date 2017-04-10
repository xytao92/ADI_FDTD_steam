//writen by liuzhichao 
//2016/10/8

#include "functions.h"
#include "GSS_ADI_matel.h"
#include "GSS_ADI_CPML.h"
#include "GSS_ADI_mur_1st.h"

int main()
{
	clock_t start_time = clock();
	initGrid(halfgrid_before, halfgrid_now);//对网格数据进行初始化
	int choice;

	while (TRUE)
	{
		cout << ">>>>>>>>>>>>>请 选 择:<<<<<<<<<<<<\n" << endl;
		
		cout << ">1. 金 属 吸 收 边 界\n" << endl;
		cout << ">2. CPML 吸 收 边 界\n" << endl;
		cout << ">3. 计 算 CFL 稳 定 性 条 件 dt 与 截 止 频 率 fc\n" << endl;
		cout << ">4. 计 算 理 论 值\n" << endl;
		cout << ">5. 结 束 程 序\n" << endl;
		cin >> choice;
		cout << '\n' << endl;
		switch (choice)
		{

		case 1:
		{
			adi_fdtd_leapforg_matel_GSS(halfgrid_now);//进行计算
			continue;
		}
		case 2:
		{
			adi_fdtd_leapforg_cpml_GSS(halfgrid_now, halfgrid_before);//进行计算
			continue;
		}
		case 3:
		{
			double result = CFL_calc();
			double fc = fc_calc();
			cout << " --------计算模型-------- \n\n" << "Z: \t" << Z << "  X: \t" << X << "  Y: \t" << Y << endl;
			cout << "\n" << endl;
			cout << " --------网格数-------- \n\n" << "Nx: \t" << Nx << "  Ny: \t" << Ny << "  Nz: \t" << Nz << endl;
			cout << "\n" << endl;
			cout << " --------网格大小-------- \n\n" << "dx: \t" << dx << "  dy: \t" << dy << "  dz: \t" << dz << endl;
			cout << "\n" << endl;
			cout << "此结构下，截止频率fc为：\t " << fc << endl;
			cout << "在满足CFL稳定性条件下: dt <= \t"<< result << endl;
			cout << "\n" << endl;
			cout << " 若要观察到至少 10个 周期的波形，输入源频率: " << " > " << 10 / (400 * dt);

			
			cout << '\n' << endl;//进行计算
			continue;
		}
		case 4:
		{
			theor_val_gen();
			cout << "理论值已产生并存储至---> "<< theor_val_filepath << "\n" << endl;
			continue;
		}
		case 5:
		{
			cout << ">>>>>>>>程序结束<<<<<<<<" << endl;
			break;
		}
		default:
		{
			cout << "抱歉，无此选项，请重新输入。\n" << endl;
			continue;
		}

		}

		free(halfgrid_before, halfgrid_now);//释放空间

		clock_t end_time = clock();
		cout << "Running time: " << static_cast<float>(end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
		system("pause");
		return 0;

	}
}