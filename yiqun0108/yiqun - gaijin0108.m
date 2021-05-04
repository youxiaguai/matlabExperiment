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
 
  MM = size(G,1);%����
    figure(3)
    axis([0,MM,0,MM])%x��y��ȡֵ��Χ
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
 %��դ���ͼתΪ�ڽӾ���
MM=size(G,1);%����G������
D=zeros(MM*MM,MM*MM);%N��ʾ����Ĺ�ģ�����ظ��������ؾ���D����
a=1;%С�������صı߳�
% Ex=a*(mod(E,MM)-0.5);%��ֹ������� a���Ա���E��MM���У�ȡ��(�õ���)���0.5 ��������
% if Ex==-0.5
%     Ex=MM-0.5;
% end
% Ey=a*(MM+0.5-ceil(E/MM));%E/MM���ȡ�� ��ֹ��������
% Eta=zeros(1,N);%����ʽ��Ϣ��ȡΪ��Ŀ����ֱ�߾���ĵ��� ��ʼ��Ϣ�ؾ���
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
%���湹������ʽ��Ϣ����
Eta=zeros(1,N);%����ʽ��Ϣ��ȡΪ��Ŀ����ֱ�߾���ĵ��� ��ʼ��Ϣ�ؾ���
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
Alpha= 1 ;                       %% Alpha������Ϣ����Ҫ�̶ȵĲ���
Beta = 6 ;                          % %  Beta��������ʽ������Ҫ�̶ȵĲ���
Rho = 0.2;  % Rho��Ϣ������ϵ��
sim= 0.3 ;     %������
Q = 1 ;                               % Q��Ϣ������ǿ��ϵ��
minkl = inf ;
minkm = inf ;
N=size(D,1);
minl = 0 ;
Tau=ones(N,N);%��Ϣ��
ROUTES=cell(K,M);%��ϸ���ṹ�洢ÿһ����ÿһֻ���ϵ�����·�� ���ϸ���*������������ÿ��Ԫ����һ���ṹ
PL=zeros(K,M);%�þ���洢ÿһ����ÿһֻ���ϵ�����·�߳���
%% -----------����K��������ʳ���ÿ���ɳ�Mֻ����--------------------
tic
for k=1:K
    %disp(k);
    for m=1:M
%%     ��һ����״̬��ʼ��
        W=S;%��ǰ�ڵ��ʼ��Ϊ��ʼ��
        Path=S;%����·�߳�ʼ��
        pathlength=0;%����·�߳��ȳ�ʼ��
        Tabu=ones(1,N); %���ɽ����б�,���нڵ��δ�߹������Զ���Ϊ1
        Tabu(S)=0;%�Ѿ��ڳ�ʼ���ˣ����Ҫ�ų�
        DD=D;%�ڽӾ����ʼ��
%%     �ڶ�������һ������ǰ���Ľڵ�
        DW=DD(W,:);
        LJD=find(DW>0);%��ѡ�ڵ㼯 �����ؿ����ߵĽڵ�����<������>
        Len_LJD=length(LJD);%���� ��ѡ�ڵ�ĸ��� 
%%     ��ʳֹͣ����������δ����ʳ�������������ͬ
        while Len_LJD>=1&&W~=N  %W~=E&&
%%         ��������ת�ֶķ�ѡ����һ����ô��
            node=zeros(1,Len_LJD); %������ѡ�ڵ�
            for i=1:Len_LJD
                 node(i)=(Tau(W,LJD(i))^Alpha)*((1/Eta(1,LJD(i)))^Beta); %w��i���ڵ�
            end
            node=node/(sum(node));%�������ʷֲ� �Ѹ���·���ĸ���ͳһ����Ϊ1��
            Pcum=cumsum(node);  %node�ۼ�ֵ 
            Select=find(Pcum>=rand);%��������0~1֮�������������̶��㷨��������������ֲ����Ž�
            to_visit=LJD(Select(1));%��һ����Ҫǰ���Ľڵ�
%%         ���Ĳ���״̬���ºͼ�¼
            Path=[Path,to_visit];%·�߽ڵ�����
            pathlength=pathlength+DD(W,to_visit);%·���������ӣ���¼���ε������·�߳��ȣ�ÿֻ���϶����Լ��߹��ĳ��ȼ�¼�������С�
            W=to_visit;%�����Ƶ���һ���ڵ�
            %N�����е�
            for n=1:N
                if Tabu(n)==0    %�����б�
                    DD(W,n)=0;  %�ڴ˴�ѭ��������Ϊ���ɴ�   
                    DD(n,n)=0;
                end
            end
            Tabu(W)=0;%�ѷ��ʹ��Ľڵ�ӽ��ɱ���ɾ��
            DW=DD(W,:);
            LJD=find(DW>0);%��ѡ�ڵ㼯
            Len_LJD=length(LJD);%��ѡ�ڵ�ĸ���
        end
%%     ���岽������ÿһ��ÿһֻ���ϵ���ʳ·�ߺ�·�߳���
        ROUTES{k,m}=Path; %��k�ε��� ��mֻ���ϵ�·��
        if Path(end)==N
            PL(k,m)=pathlength; %����Ŀ����·�߳���
        else
            PL(k,m)=inf;  %��������ͬ
        end
    end
%% ��������������Ϣ��
    Delta_Tau=zeros(N,N);%��������ʼ��
    for mm=1:M %Mֻ����
        if PL(k,mm)<inf %˳������Ŀ��������·�߳���
            ROUT=ROUTES{k,mm}; %����·��
            TS=length(ROUT)-1; %���� ����ת�ƴ���
            PL_km=PL(k,mm);%·�߳���
            for s=1:TS
                x=ROUT(s); %��һ���ڵ�
                y=ROUT(s+1); %��һ���ڵ�
                Delta_Tau(x,y)=Delta_Tau(x,y)+Q/PL_km; %(x,y)�������ڵ�֮��Ĺ�ϵ(��Ϣ����) ϵ������·�߳���
                Delta_Tau(y,x)=Delta_Tau(y,x)+Q/PL_km;
            end
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau;%��Ϣ�ػӷ�һ���֣�������һ����
end
toc
%% ---------------------------��ͼ--------------------------------
figure(4)
axis([0,K,0,20])
for i=1:K
   ggb=PL(i,:);
   mmmda(i)=sum(ggb)/50;
end
plot(mmmda)
plotif=1; %�Ƿ��ͼ�Ŀ��Ʋ���
if plotif==1
    %����������
   
    meanPL=zeros(1,K); %k����������
    minPL=zeros(1,K);
    for i=1:K
        PLK=PL(i,:); %����i�ε�������·�߳��ȸ�ֵ��PLK
        Nonzero=find(PLK<inf);%����һϵ�п���·�ߵı��
        if length(Nonzero)~=0
            PLKPLK=PLK(Nonzero);%���¿���·�ߣ���������
            meanPL(i)=mean(PLKPLK); %��ȡ��ο���·����ƽ��ֵ
            minPL(i)=min(PLKPLK);%�����С·��
        end
    end
end

plotif3=1;%�������������ͼ
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
        PLK=PL(k,:); %����k�ε�������·�߳��ȸ�ֵ��PLK
        minPLK=min(PLK);
        if(minPLK<minmumPLK)
        pos=find(PLK==minPLK); %�ҵ����������·�߳�����ȵ�·�����
        minmumPLK=minPLK;
        minm=pos(1);
        mink=k ; %����k��
        end
    end
        ROUT=ROUTES{mink,minm}; %�ҳ���С·��·��
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



