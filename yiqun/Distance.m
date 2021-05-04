function D = Distance(citys)
%% ������������֮��ľ���
% ���룺�����е�λ�����꣨citys��
% �������������֮��ľ��루D��

n = size(citys,1);
D = zeros(n,n);
for i = 1:n
    for j = i + 1:n
        D(i,j)=sqrt((citys(i,1)-citys(j,1))^2 + (citys(i,2) - citys(j,2))^2);
        D(j,i)=D(i,j);
    end
    D(i,i)=1e-4;        %�Խ��ߵ�ֵΪ0�������ں������������Ҫȡ�����������һ����С������0
end