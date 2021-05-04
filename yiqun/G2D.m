function D=G2D(G)   %相邻栅格点初始化
l=size(G,1);         %l=G矩阵的行数
%D是所有点到所有点的矩阵。纵向将所有点作为起点，横向将所有点作为邻域局部到达点
D=zeros(l*l,l*l); 
for i=1:l 
    for j=1:l 
        if G(i,j)==0 %判断局部起始栅格是否有障碍
            for m=1:l 
                for n=1:l 
                    if G(m,n)==0 %判断局部终点栅格，即相邻栅格是否有障碍
                        im=abs(i-m);jn=abs(j-n); 
                        if im+jn==1||(im==1&&jn==1) %判断该栅格是否为当前栅格的邻域栅格
                        D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5; %将所有点到所有邻域点的代价值填入D中
                                          %D((i-1)*l+j,(m-1)*l+n)在矩阵中定位邻域矩阵的具体位置
                                          %(im+jn)^0.5在具体位置中填入具体代价值 
                        end 
                    end 
                end 
            end 
        end 
    end 
end
