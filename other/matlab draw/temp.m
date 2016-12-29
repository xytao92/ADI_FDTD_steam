  A=load('E:\ADI_FDTD\ADI_FDTD_steam\result\temp_mur_2nd.txt'); %将sj.txt文件数据，赋值给A矩阵
   %A=load('E:\ADI_FDTD\Source-test\source_test\result\source.txt'); %将sj.txt文件数据，赋值给A矩阵
 x1=A(:,1); 
 y1=A(:,2);
 y2=A(:,3);
 y3=A(:,4);

plot(x1,y1,'LineWidth',2,'Color',[1 0 0],'DisplayName','Ex');
hold on
plot(x1,y2,'LineWidth',2,'Color',[0.0784313753247261 0.168627455830574 0.549019634723663],'DisplayName','Ey');
hold on
plot(x1,y3,'LineWidth',2,'Color',[0 1 0],'DisplayName','Ez');
 
 % 创建 ylabel
ylabel('输出');

% 创建 xlabel
xlabel('时间步长(1.0*e-12s)');

% 创建 title
title('小网格（50*20*10）下输出结果','FontSize',12,'FontName','Microsoft YaHei UI');
% 创建 legend
legend('show');
grid on

