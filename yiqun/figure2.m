figure(2) %建立图形窗口
   axis([0,MM,0,MM]) %设置图的横纵坐标，MM为地图矩阵的行数或列数
   for i=1:MM 
      for j=1:MM 
       %此处的行纵坐标从0开始，为对应于矩阵地图从1开始且pcolor加了一行，所以此处横纵坐标要换算
       %具体操作之后在优化
         if G(i,j)==1           %0是黑色代表障碍，1为白色无障碍
            x1=j-1;y1=MM-i; 
            x2=j;y2=MM-i; 
            x3=j;y3=MM-i+1; 
            x4=j-1;y4=MM-i+1; 
            fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]); %将1234点所围成的图形进行白色填充，fill为颜色填充
            hold on 
         else 
            x1=j-1;y1=MM-i; 
            x2=j;y2=MM-i; 
            x3=j;y3=MM-i+1; 
            x4=j-1;y4=MM-i+1; 
            fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]); 
            hold on 
           set(gca,'YDir','reverse');%地图的矩阵的行列与正常坐标轴的行列一致
                                   % set(gca,'propertyname','propertyvalue'......)命令可以调整图形的坐标属性
                                   % ydir：y轴方向；reverse：反转
         end
      end
   end
   hold on
   title('机器人运动轨迹'); 
   xlabel('坐标x'); 
   ylabel('坐标y');
   ROUT=ROUTES{mink,minl}; %ROUT取最短行走路线，mink为该路线的迭代号,minl为该路线的蚂蚁号
   LENROUT=length(ROUT); %提取路线长度，即走了几个栅格
   Rx=ROUT; %Rx与Ry中分别存储具体的该路线
   Ry=ROUT; 
%将该路线的栅格索引号转换为横纵坐标
   for ii=1:LENROUT 
     Rx(ii)=a*(mod(ROUT(ii),MM)-0.5); 
     if Rx(ii)==-0.5 
        Rx(ii)=MM-0.5; 
     end 
     Ry(ii)=a*(MM+0.5-ceil(ROUT(ii)/MM)); 
   end 
plot(Rx,Ry) %绘制路线
