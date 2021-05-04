%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�������ƣ�ѡ���� selection.m
%%˵����
function [newx,newy]=selection(x,y,fitvalue)

totalfit=sum(fitvalue); %����Ӧֵ֮��
fitvalue1=fitvalue/totalfit; %�������屻ѡ��ĸ���
fitvalue=cumsum(fitvalue1); %�� fitvalue=[1 2 3 4]���� cumsum(fitvalue)=[1 3 6 10] 
[px,py]=size(fitvalue);
ms=sort(rand(px,1)); %��С��������       ����һ�� px�� 1�� ���������Ȼ���С��������
fitin=1;
newin=1;
while newin<=px
if (ms(newin))<fitvalue(fitin) && fitvalue1(fitin)>0  %fitvalue1(fitin)>0 ��֤��Ӧ��Ϊ0�ĸ��岻��ѡ��
newx(newin,:)=x(fitin,:);
newy(newin,:)=y(fitin,:);
newin=newin+1;
else
fitin=fitin+1;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�������ƣ����и���ѡ���� best.m
%%��ڲ�������ʼ��Ⱥ    ��Ⱥ��Ӧ��
%%���ڲ�������Ѹ���  ��Ѹ�����Ӧ��
%%˵����
    %%������Ӧ�ȴ�С����ѡ����Ѹ��塣��Ӧ�����ĸ��彫��ѡ������Ϊ��������ֵ���ء�
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
