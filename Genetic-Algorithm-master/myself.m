clear;clc;close all;

%%�Ŵ���������
NUMPOP=100;%��ʼ��Ⱥ��С
irange_l=-1; %���������
irange_r=2;
LENGTH=22; %�����Ʊ��볤��
ITERATION = 10000;%��������
CROSSOVERRATE = 0.7;%�ӽ���
SELECTRATE = 0.5;%ѡ����
VARIATIONRATE = 0.001;%������

%��ʼ����Ⱥ
pop=m_InitPop(NUMPOP,irange_l,irange_r);
pop_save=pop;
%���Ƴ�ʼ��Ⱥ�ֲ�
x=linspace(-1,2,1000);
y=m_Fx(x);
plot(x,y);
hold on
for i=1:size(pop,2)
    plot(pop(i),m_Fx(pop(i)),'ro');
end
hold off
title('��ʼ��Ⱥ');

%��ʼ����
for time=1:ITERATION
    %�����ʼ��Ⱥ����Ӧ��
    fitness=m_Fitness(pop);
    %ѡ��
    pop=m_Select(fitness,pop,SELECTRATE);
    %����
    binpop=m_Coding(pop,LENGTH,irange_l);
    %����
    kidsPop = crossover(binpop,NUMPOP,CROSSOVERRATE);
    %����
    kidsPop = Variation(kidsPop,VARIATIONRATE);
    %����
    kidsPop=m_Incoding(kidsPop,irange_l);
    %������Ⱥ
    pop=[pop kidsPop];
end
figure
x=linspace(-1,2,1000);
y=m_Fx(x);
plot(x,y);
hold on
for i=1:size(pop,2)
    plot(pop(i),m_Fx(pop(i)),'ro');
end
hold off
title('��ֹ��Ⱥ');

disp(['���Ž⣺' num2str(max(m_Fx(pop)))]);
disp(['�����Ӧ�ȣ�' num2str(max(m_Fitness(pop)))]);   
    
    


