%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%函数名称：选择函数 selection.m
%%说明：
function [newx,newy]=selection(x,y,fitvalue)

totalfit=sum(fitvalue); %求适应值之和
fitvalue1=fitvalue/totalfit; %单个个体被选择的概率
fitvalue=cumsum(fitvalue1); %如 fitvalue=[1 2 3 4]，则 cumsum(fitvalue)=[1 3 6 10] 
[px,py]=size(fitvalue);
ms=sort(rand(px,1)); %从小到大排列       生成一个 px行 1列 的随机矩阵，然后从小到大排列
fitin=1;
newin=1;
while newin<=px
if (ms(newin))<fitvalue(fitin) && fitvalue1(fitin)>0  %fitvalue1(fitin)>0 保证适应度为0的个体不被选中
newx(newin,:)=x(fitin,:);
newy(newin,:)=y(fitin,:);
newin=newin+1;
else
fitin=fitin+1;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%函数名称：最有个体选择函数 best.m
%%入口参数：初始种群    种群适应度
%%出口参数：最佳个体  最佳个体适应度
%%说明：
    %%按照适应度大小进行选择最佳个体。适应度最大的个体将被选出，作为函数返回值返回。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bestx,besty,bestfit]=best(x,y,fitvalue)
    [px,py]=size(x);
    bestx=x(1,:);
    besty=y(1,:);
    bestfit=fitvalue(1);
    for i=2:px
        if fitvalue(i)>bestfit
            bestx=x(i,:);
            besty=y(i,:);
            bestfit=fitvalue(i);
        end
    end
end
