function fitness=m_Fitness(pop)
%% Fitness Function
%y=xsin(3x)��[-1,2]�ϣ����ֵҲ���ᳬ��2
%���Լ��㺯��ֵ��2�ľ��룬������Сʱ����Ϊ���Ž�
%��Ӧ�Ⱥ���Ϊ1/����
for n=1:size(pop,2)
    fitness(n)=1/(2-m_Fx(pop(:,n)));
end

end
