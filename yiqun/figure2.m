figure(2) %����ͼ�δ���
   axis([0,MM,0,MM]) %����ͼ�ĺ������꣬MMΪ��ͼ���������������
   for i=1:MM 
      for j=1:MM 
       %�˴������������0��ʼ��Ϊ��Ӧ�ھ����ͼ��1��ʼ��pcolor����һ�У����Դ˴���������Ҫ����
       %�������֮�����Ż�
         if G(i,j)==1           %0�Ǻ�ɫ�����ϰ���1Ϊ��ɫ���ϰ�
            x1=j-1;y1=MM-i; 
            x2=j;y2=MM-i; 
            x3=j;y3=MM-i+1; 
            x4=j-1;y4=MM-i+1; 
            fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]); %��1234����Χ�ɵ�ͼ�ν��а�ɫ��䣬fillΪ��ɫ���
            hold on 
         else 
            x1=j-1;y1=MM-i; 
            x2=j;y2=MM-i; 
            x3=j;y3=MM-i+1; 
            x4=j-1;y4=MM-i+1; 
            fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]); 
            hold on 
           set(gca,'YDir','reverse');%��ͼ�ľ�������������������������һ��
                                   % set(gca,'propertyname','propertyvalue'......)������Ե���ͼ�ε���������
                                   % ydir��y�᷽��reverse����ת
         end
      end
   end
   hold on
   title('�������˶��켣'); 
   xlabel('����x'); 
   ylabel('����y');
   ROUT=ROUTES{mink,minl}; %ROUTȡ�������·�ߣ�minkΪ��·�ߵĵ�����,minlΪ��·�ߵ����Ϻ�
   LENROUT=length(ROUT); %��ȡ·�߳��ȣ������˼���դ��
   Rx=ROUT; %Rx��Ry�зֱ�洢����ĸ�·��
   Ry=ROUT; 
%����·�ߵ�դ��������ת��Ϊ��������
   for ii=1:LENROUT 
     Rx(ii)=a*(mod(ROUT(ii),MM)-0.5); 
     if Rx(ii)==-0.5 
        Rx(ii)=MM-0.5; 
     end 
     Ry(ii)=a*(MM+0.5-ceil(ROUT(ii)/MM)); 
   end 
plot(Rx,Ry) %����·��
