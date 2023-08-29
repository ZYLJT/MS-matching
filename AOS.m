function [U]=AOS(U_prev,step,diff_c)
%�ú���ʵ�ּ��Է����㷨
%U_prev��ǰһ��ĸ�������ͼ��
%step�����õĲ������ӣ��������ĳ߶���sigma(i),��ǰһ��ĳ߶���sigma(i-1)
%��˲���������1/2*(sigma(i)^2-sigma(i-1)^2)
%diff_c�Ǹ���ǰһ��߶�ͼ���˹�˲������õ�����ɢϵ����СҲ��M*N
%U�ǵ�ǰ��ĸ������Գ߶�ͼ�񣬴�С��U_prevһ��

[U1]=AOS_row(U_prev,step,diff_c);
[U2]=AOS_col(U_prev,step,diff_c);
% U=1/2*(U1+U2);
U=1/2*(U1+U2);   %��ʾ��L=L1+L2

end


function [U1]=AOS_row(U1_prev,step,diff_c)
%�ú������з��������ɢ
%%
[M,N]=size(U1_prev);
U1=zeros(M,N);
%����ÿһ��
for i=1:1:M
    d=U1_prev(i,:);%������
    
    %���Ǿ���ĶԽ��߲���
    a=diff_c(i,:);
    a(2:N-1)=2*a(2:N-1);
    a(1:N-1)=a(1:N-1)+diff_c(i,2:N);
    a(2:N)=a(2:N)+diff_c(i,1:N-1);
    a=-1/2*a;
    
    %���Ǿ�������
    b=diff_c(i,1:N-1)+diff_c(i,2:N);
    b=1/2*b;
    
    %���Ǿ�������,�þ���Գƣ����c=b
    c=b;
    
    %�������ԽǷ�����Ľ�
    a=1-2*step*a;
    b=-2*step*b;
    c=-2*step*c;
    x=thomas_algorith(a,b,c,d);   %��ҪĿ���ǽ��������A
    U1(i,:)=x;
end

end
    
function [U2]=AOS_col(U2_prev,step,diff_c)
%�ú������з������ɢ
%%
[M,N]=size(U2_prev);
U2=zeros(M,N);
%����ÿһ��
for i=1:1:N
    d=U2_prev(:,i);%������
    %���Ǿ���ĶԽ��߲���
    a=diff_c(:,i);
    a(2:M-1)=2*a(2:M-1);
    a(1:M-1)=a(1:M-1)+diff_c(2:M,i);
    a(2:M)=a(2:M)+diff_c(1:M-1,i);
    a=-1/2*a;
    
    %���Ǿ�������
    b=diff_c(1:M-1,i)+diff_c(2:M,i);
    b=1/2*b;
    
    %���Ǿ�������,�þ���Գƣ����c=b
    c=b;
    
    %�������ԽǷ�����Ľ�
    a=1-2*step*a;
    b=-2*step*b;
    c=-2*step*c;
    x=thomas_algorith(a',b',c',d');
    U2(:,i)=x';
end
end

%% 
function [x]=thomas_algorith(a,b,c,d)
%�ú�������thomas�����Է�����Ľ⣬�����������������Խ�����ϵͳ
%���������ʽ��,ע��a,b,��r��λ��
%[a1,b1,0,0..............0,0,0][x1]   =  [d1]
%[c1,a2,b2,0,0...........0,0,0][x2]   = [d2]
%[0,c2,a3,b3,0,..........0,0,0][x3]    =[d3]
%[0,.....................,,,,,][x.]    = [d.]
%[0,0,0,0,0,....cM-2,aM-1,bM-1][xM-1]  = [dM-1]
%[0,0,0,0,0,....0,0,0,cM-1,aM ][xM]    = [dM]
%����a��d��С��Mά��������b��c��С��M-1ά������

%% LR�ֽ⣬��LR=A
%L=[1  0  0  0  0 . . . 0]
%  [L1 1  0  0  0.......0]
%  [0  L2 0  0  0 ......0]
%  [0  0   LM-3  1   0  0]
%  [0  0        LM-2 1  0]
%  [0  0           LM-1 1]

%R=[m1  r1  0  0  0 ......0]
%  [0   m2  r2 0  0.......0]
%  [0   0   m3 r3 0.......0]
%  [0   0    mM-2  rM-2...0]
%  [0   0        mM-1  rM-1]
%  [0   0          0     mM]

[~,N]=size(a);%����b��һ��������
m=zeros(1,N);
L=zeros(1,N-1);

%% ����r����֪r=b;
% r=b;

%% ����m��L
m(1)=a(1);
for i=1:1:N-1
    L(i)=c(i)/m(i);
    m(i+1)=a(i+1)-L(i)*b(i);
end

%% LRx=d,��Ly=d,����y
y=zeros(1,N);
y(1)=d(1);
for i=2:1:N
    y(i)=d(i)-L(i-1)*y(i-1);
end

%% Rx=y
x=zeros(1,N);
x(N)=y(N)/m(N);
for i=N-1:-1:1
    x(i)=(y(i)-b(i)*x(i+1))/m(i);
end

end

