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
%%�������ƣ���ʼ����Ⱥ popinit.m
%%��ڲ�������Ⱥ����  ��������
%%���ڲ�������ʼ��Ⱥ
%%˵����
    %%��ʼ��Ⱥ�ĸ����X��������Y������ֿ���ţ��ֱ���ھ��� x,y�У���Ϊ��������ֵ����
    %%��ʼ��Ⱥ�Ĳ�������ȥ��ʼ������ֹ�����㣬�������x�ᡢy��������������Ӵ�С��������
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
