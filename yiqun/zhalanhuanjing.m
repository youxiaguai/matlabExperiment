%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%建立环境地图%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zhalanhuanjing(map)
n = size(map);
step = 1;
a = 0 : step :n(1);
b = 0 : step :n(2);
figure(1)
axis([0 n(2) 0 n(1)]); %设置地图横纵尺寸
set(gca,'xtick',b,'ytick',a,'GridLineStyle','-',...
'xGrid','on','yGrid','on');
hold on
r = 1;
for(i=1:n(1))         %设置障碍物的左下角点的x,y坐标
    for(j=1:n(2))
        if(map(i,j)==1)
            p(r,1)=j-1;
            p(r,2)=i-1;
            fill([p(r,1) p(r,1) + step p(r,1) + step p(r,1)],...
                 [p(r,2) p(r,2) p(r,2) + step p(r,2) + step ],'k');
            r=r+1;
            hold on
        end
    end
end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%栅格数字标识%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_text = 1:1:n(1)*n(2); %产生所需数值.
for i = 1:1:n(1)*n(2)
    [row,col] = ind2sub([n(2),n(1)],i);
    text(row-0.9,col-0.5,num2str(x_text(i)),'FontSize',8,'Color','0.7 0.7 0.7');
end
hold on
axis square
