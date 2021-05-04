clear all
clc
G=[ 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 1 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 1 1 0 0 0 1 0 1 1 0 1 0 0 1 0 1 0 0 
    1 0 0 1 0 0 1 0 1 0 0 0 1 0 1 0 0 1 0 0 
    0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 1 0 0 0 0 
    0 1 0 1 1 1 0 0 1 1 0 0 0 1 0 0 1 0 0 0 
    0 1 0 1 0 0 0 0 0 1 1 1 1 0 0 0 0 0 1 0 
    0 0 0 0 1 0 1 0 0 0 1 0 0 0 0 1 1 1 0 0 
    0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 
    0 0 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 0
    0 0 0 0 0 1 0 1 1 1 0 0 0 1 0 0 0 0 0 0
    0 0 0 0 1 0 1 0 0 0 0 0 1 0 1 0 0 0 0 0
    0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 1 1 0 0 0
    0 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 1 1 0
    1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 0 0 0 1
    0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 1 0 0
    0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 0 0
    0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0
    ];
 
  MM = size(G,1);%行数
    figure(3)
    axis([0,MM,0,MM])%x，y的取值范围
    for i=1:MM
        for j=1:MM
            if G(i,j)==1
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.3,0.3,0.3]);
                hold on
            else
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]);
                hold on
           end
        end
    end
 %把栅格地图转为邻接矩阵
MM=size(G,1);%返回G的行数
D=zeros(MM*MM,MM*MM);%N表示问题的规模（象素个数）返回矩阵D行数
a=1;%小方格象素的边长
% Ex=a*(mod(E,MM)-0.5);%终止点横坐标 a乘以变量E对MM（行）取余(得到列)后减0.5 即所处列
% if Ex==-0.5
%     Ex=MM-0.5;
% end
% Ey=a*(MM+0.5-ceil(E/MM));%E/MM结果取整 终止点纵坐标
% Eta=zeros(1,N);%启发式信息，取为至目标点的直线距离的倒数 初始信息素矩阵
for i= 1:MM
   for j=1:MM 
       if G(i,j)==0
           for m=1:MM
               for n=1:MM
                   if G(m,n)==0 && ((i-m)^2 + (j-n)^2)<=2
                       im=abs(i-m);
                       jn=abs(j-n);
                       if(im+jn==1)||(im==1&&jn==1)
                           D((i-1)*MM+j,(m-1)*MM+n)=(im+jn)^0.5;
                       end                       
                   end
               end
           end
       end
   end
end
for i= 1:MM
   for j=1:MM 
       if G(i,j)==1
           for m=1:MM
               for n=1:MM
                   if G(m,n)==1
                       im=abs(i-m);
                       jn=abs(j-n);
                       if(im==1&&jn==1)
                           if j>n
                               if i>m
                               D((i-2)*MM+j,(i-1)*MM+j-1)=0;  %YY=(i-2)*MM+j;
                               end
                               if i<m
                                D(i*MM+j,(i-1)*MM+j-1)=0; %YY=i*MM+j;
                               end
                                     %%XX=(i-1)*MM+j-1;
                           end
                           if j<n
                               %XX=(i-1)*MM+j+1;
                               if i>m
                               D((i-2)*MM+j,(i-1)*MM+j+1)=0;  %YY=(i-2)*MM+j;
                               end
                               if i<m
                                D(i*MM+j,(i-1)*MM+j+1)=0; %YY=i*MM+j;
                               end
                           end
                           
                       end                       
                   end
               end
           end
       end
   end
end
N=size(D,1);
%下面构造启发式信息矩阵
Eta=zeros(1,N);%启发式信息，取为至目标点的直线距离的倒数 初始信息素矩阵
for i=1:MM
    for j=1:MM
        Eta(1,(i-1)*MM+j)=((MM-i)^2+(MM-j)^2)^0.5;      
    end
end
Eta(1,MM*MM)=0.01;
K=200;
S=1;
M=100;
p = 2 ;
Alpha= 1 ;                       %% Alpha表征信息素重要程度的参数
Beta = 6 ;                          % %  Beta表征启发式因子重要程度的参数
Rho = 0.2;  % Rho信息素蒸发系数
sim= 0.3 ;     %西格玛
Q = 1 ;                               % Q信息素增加强度系数
minkl = inf ;
minkm = inf ;
N=size(D,1);
minl = 0 ;
Tau=ones(N,N);%信息素
ROUTES=cell(K,M);%用细胞结构存储每一代的每一只蚂蚁的爬行路线 蚂蚁个数*迭代次数矩阵，每个元素是一个结构
PL=zeros(K,M);%用矩阵存储每一代的每一只蚂蚁的爬行路线长度
%% -----------启动K轮蚂蚁觅食活动，每轮派出M只蚂蚁--------------------
tic
for k=1:K
    %disp(k);
    for m=1:M
%%     第一步：状态初始化
        W=S;%当前节点初始化为起始点
        Path=S;%爬行路线初始化
        pathlength=0;%爬行路线长度初始化
        Tabu=ones(1,N); %生成禁忌列表,所有节点均未走过，所以都置为1
        Tabu(S)=0;%已经在初始点了，因此要排除
        DD=D;%邻接矩阵初始化
%%     第二步：下一步可以前往的节点
        DW=DD(W,:);
        LJD=find(DW>0);%可选节点集 即返回可以走的节点坐标<矩阵编号>
        Len_LJD=length(LJD);%计数 可选节点的个数 
%%     觅食停止条件：蚂蚁未遇到食物或者陷入死胡同
        while Len_LJD>=1&&W~=N  %W~=E&&
%%         第三步：转轮赌法选择下一步怎么走
            node=zeros(1,Len_LJD); %遍历可选节点
            for i=1:Len_LJD
                 node(i)=(Tau(W,LJD(i))^Alpha)*((1/Eta(1,LJD(i)))^Beta); %w行i个节点
            end
            node=node/(sum(node));%建立概率分布 把各个路径的概率统一到和为1；
            Pcum=cumsum(node);  %node累计值 
            Select=find(Pcum>=rand);%产生任意0~1之间的随机数，轮盘赌算法，尽量避免陷入局部最优解
            to_visit=LJD(Select(1));%下一步将要前往的节点
%%         第四步：状态更新和记录
            Path=[Path,to_visit];%路线节点增加
            pathlength=pathlength+DD(W,to_visit);%路径长度增加，记录本次迭代最佳路线长度，每只蚂蚁都有自己走过的长度记录在向量中。
            W=to_visit;%蚂蚁移到下一个节点
            %N：所有点
            for n=1:N
                if Tabu(n)==0    %禁忌列表
                    DD(W,n)=0;  %在此次循环中设置为不可达   
                    DD(n,n)=0;
                end
            end
            Tabu(W)=0;%已访问过的节点从禁忌表中删除
            DW=DD(W,:);
            LJD=find(DW>0);%可选节点集
            Len_LJD=length(LJD);%可选节点的个数
        end
%%     第五步：记下每一代每一只蚂蚁的觅食路线和路线长度
        ROUTES{k,m}=Path; %第k次迭代 第m只蚂蚁的路线
        if Path(end)==N
            PL(k,m)=pathlength; %到达目标点的路线长度
        else
            PL(k,m)=inf;  %进入死胡同
        end
    end
%% 第六步：更新信息素
    Delta_Tau=zeros(N,N);%更新量初始化
    for mm=1:M %M只蚂蚁
        if PL(k,mm)<inf %顺利到达目标点的蚂蚁路线长度
            ROUT=ROUTES{k,mm}; %具体路线
            TS=length(ROUT)-1; %跳数 蚂蚁转移次数
            PL_km=PL(k,mm);%路线长度
            for s=1:TS
                x=ROUT(s); %上一个节点
                y=ROUT(s+1); %下一个节点
                Delta_Tau(x,y)=Delta_Tau(x,y)+Q/PL_km; %(x,y)即两个节点之间的关系(信息素量) 系数除以路线长度
                Delta_Tau(y,x)=Delta_Tau(y,x)+Q/PL_km;
            end
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau;%信息素挥发一部分，新增加一部分
end
toc
%% ---------------------------绘图--------------------------------
figure(4)
axis([0,K,0,20])
for i=1:K
   ggb=PL(i,:);
   mmmda(i)=sum(ggb)/50;
end
plot(mmmda)
plotif=1; %是否绘图的控制参数
if plotif==1
    %绘收敛曲线
   
    meanPL=zeros(1,K); %k：迭代次数
    minPL=zeros(1,K);
    for i=1:K
        PLK=PL(i,:); %将第i次迭代爬行路线长度赋值给PLK
        Nonzero=find(PLK<inf);%返回一系列可行路线的编号
        if length(Nonzero)~=0
            PLKPLK=PLK(Nonzero);%留下可行路线，重新排列
            meanPL(i)=mean(PLKPLK); %求取这次可行路径的平均值
            minPL(i)=min(PLKPLK);%提出最小路径
        end
    end
end

plotif3=1;%绘最短蚂蚁爬行图
if plotif3==1
    figure(2)
    axis([0,MM,0,MM])
    for i=1:MM
        for j=1:MM
            if G(i,j)==1
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]);
                hold on
            else
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]);
                hold on
            end
        end
    end
    minmumPLK=inf;
    for k=1:K
        PLK=PL(k,:); %将第k次迭代爬行路线长度赋值给PLK
        minPLK=min(PLK);
        if(minPLK<minmumPLK)
        pos=find(PLK==minPLK); %找到与最短爬行路线长度相等的路径标号
        minmumPLK=minPLK;
        minm=pos(1);
        mink=k ; %迭代k次
        end
    end
        ROUT=ROUTES{mink,minm}; %找出最小路径路线
        LENROUT=length(ROUT);
        Rx=ROUT;
        Ry=ROUT;
    for ii=1:LENROUT
            Rx(ii)=a*(mod(ROUT(ii),MM)-0.5);
            if Rx(ii)==-0.5
                Rx(ii)=MM-0.5;
            end
            Ry(ii)=a*(MM+0.5-ceil(ROUT(ii)/MM));
     end
        plot(Rx,Ry)
        hold on
end



