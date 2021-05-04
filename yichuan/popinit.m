function [newx,newy]=mutation(x,y,pm)

[px,py]=size(x);
newx=x;
newy=y;
for i=1:1:px
if(rand<pm)
mpoint=round(rand*py);

if mpoint<=1
    mpoint=2;
end
if mpoint==py
    mpoint=py-1;
end

newx(i,mpoint)=round(rand*py);
newy(i,mpoint)=round(rand*py);
end
end
for i=1:1:px
newx(i,:)=sort(newx(i,:));
newy(i,:)=sort(newy(i,:));
end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%函数名称：初始化族群 popinit.m
%%入口参数：种群数量  基因数量
%%出口参数：初始种群
%%说明：
    %%初始种群的个点的X轴坐标与Y轴坐标分开存放，分别放在矩阵 x,y中，作为函数返回值返回
    %%初始种群的产生，除去起始点与终止点两点，其他点的x轴、y轴随机产生，并从大到小进行排列
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y]=popinit(popsize,chromlength)

x=20.0*rand(popsize,chromlength);
y=20.0*rand(popsize,chromlength);
x(:,1)=0;
y(:,1)=0;
x(:,chromlength)=20;
y(:,chromlength)=20;
[px,py]=size(x);
for i=1:1:px
x(i,:)=sort(x(i,:));
y(i,:)=sort(y(i,:));
end 
