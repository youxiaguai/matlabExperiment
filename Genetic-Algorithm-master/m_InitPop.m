function pop=m_InitPop(numpop,irange_l,irange_r)
%% ��ʼ����Ⱥ
%  ���룺numpop--��Ⱥ��С��
%       [irange_l,irange_r]--��ʼ��Ⱥ���ڵ�����
pop=[];
for i=1:numpop
    pop(:,i)=irange_l+(irange_r-irange_l)*rand;
end
end
    