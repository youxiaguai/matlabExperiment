clear
clc
close all

f = @(x) sin(x) + x .* cos(x);   % �������ʽ
ezplot(f, [0, 2*pi])             % ��������ͼ��

N = 50;                          % ��Ⱥ����
ger = 100;                       % ��������
L = 5;                           % ���򳤶�
pc = 0.8;                        % �������
pm = 0.1;                        % �������
dco = [10000; 1000; 100; 10 ;1]; % ������
dna = randi([0, 9], [N, L]);     % ����
hold on
x = dna * dco / 99999 * 2 * pi;  % �Գ�ʼ��Ⱥ����
plot(x, f(x),'ko','linewidth',3) % ������ʼ���λ��

x1 = zeros(N, L);                % ��ʼ���Ӵ�����������
x2 = x1;                         % ͬ��
x3 = x1;                         % ͬ��
fi = zeros(N, 1);                % ��ʼ����Ӧ�ȣ�����

for epoch = 1: ger               % ��������Ϊ100
    for i = 1: N                 % �������
        if rand < pc
           d = randi(N);            % ȷ����һ������ĸ���
           m = dna(d,:);            % ȷ����һ������ĸ���
           d = randi(L-1);          % ȷ������ϵ�
           x1(i,:) = [dna(i,1:d), m(d+1:L)];  % �¸��� 1        
           x2(i,:) = [m(1:d), dna(i,d+1:L)];  % �¸��� 2
        end
    end
    x3 = dna;
    for i = 1: N                           % �������
        if rand < pm
            x3(i,randi(L)) = randi([0, 9]);
        end
    end
    dna = [dna; x1; x2; x3];               % �ϲ��¾ɻ���
    fi = f(dna * dco / 99999 * 2 * pi);    % ������Ӧ�ȣ��������
    dna = [dna, fi];
    dna = flipud(sortrows(dna, L + 1));    % ����Ӧ�Ƚ�������
    while size(dna, 1) > N                 % ��Ȼѡ��
        d = randi(size(dna, 1));           % ������
        if rand < (d - 1) / size(dna, 1)
            dna(d,:) = [];
            fi(d, :) = [];
        end
    end
    dna = dna(:, 1:L);
end
x = dna * dco / 99999 * 2 * pi;            % ��������Ⱥ����
plot(x, f(x),'ro','linewidth',3)           % �������ս��λ��
disp(['���Ž�Ϊx=',num2str(x(1))]);
disp(['����ֵΪy=',num2str(fi(1))]);