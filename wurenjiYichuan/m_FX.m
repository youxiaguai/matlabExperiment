function J=m_FX(Jthreat,Jfuel)
%% 要求解的函数
    J=u*Jthreat+(1-u)*Jfuel;
end