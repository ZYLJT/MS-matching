%  create_nonelinear_space:�������������Գ߶ȿռ�
%  Input:  
%           im��       �������ԭʼͼ��;
%           layers���ǹ����ĳ߶ȿռ�Ĳ���������û��ʹ���²�������.
%           scale_value��Ӱ��ĳ߶����ȱ���ϵ����Ĭ����1.6;
%           which_diff�� ������ʹ���ĸ��������㴫��ϵ��ȡֵ��1,2,3;
%           sigma_1���ǵ�һ���ͼ��ĳ߶ȣ�Ĭ����1.6���߶ȿռ��һ���ͼ����image������׼��
%           sigma_2����ÿ�μ�����һ��ͼ��֮ǰ����֮ǰ��ͼ��ĸ�˹ƽ����׼��,Ĭ����1����
%           ratio�����������ĳ߶ȱ�
%           perc���Ǽ���Աȶ����ӵ��ݶȰٷ�λ������Ĭ����0.7
%  Output:
%              Nonelinear_Scalespace���ǹ����ķ�����Ӱ��߶ȿռ�

function [Nonelinear_Scalespace]=HAPCG_nonelinear_space(im,layers,scale_value,which_diff,...
                                                        sigma_1,sigma_2,...
                                                        ratio,...
                                                        perc )
%% Ĭ�ϲ�������
if nargin < 2
    layers               = 3;          %  ������Ӱ�����������.  
end
if nargin < 3
    scale_value       = 1.6;       %  Ӱ��߶�����ϵ��   
end
if nargin < 4
    which_diff         = 2;         %  ������ʹ���ĸ��������㴫��ϵ��ȡֵ��1,2,3   
end
if nargin < 5
    sigma_1           = 1.6;         %  ��һ���ͼ��ĳ߶ȣ�Ĭ����1.6.
end
if nargin < 6
    sigma_2           = 1;            %  ��ÿ�μ�����һ��ͼ��֮ǰ����֮ǰ��ͼ��ĸ�˹ƽ����׼��.
end
if nargin < 7
     ratio               = 2^(1/3);   %  ��������ĳ߶ȱ�
end
if nargin < 8
     perc               = 0.7;            %  ����Աȶ����ӵ��ݶȰٷ�λ.
end
%% ��Ӱ��ת��Ϊ�Ҷ�ͼ
[~,~,num1]=size(im);
if(num1==3)
    dst=rgb2gray(im);
else
    dst=im;
end
% ��Ӱ��ת��Ϊ������Ӱ����ֵ��[0~1]֮�� 
image=im2double(dst);
[M,N]=size(image);
%% ��ʼ��������Ӱ��cell�ռ�
Nonelinear_Scalespace=cell(1,layers);
for i=1:1:layers
    Nonelinear_Scalespace{i}=zeros(M,N);
end

%���ȶ�����ͼ����и�˹ƽ��
windows_size=2*round(2*sigma_1)+1;
W=fspecial('gaussian',[windows_size windows_size],sigma_1);      % Fspecial�������ڴ���Ԥ������˲�����
image=imfilter(image,W,'replicate');                                              %base_image�ĳ߶���sigma_1  % ����������������άͼ������˲���
Nonelinear_Scalespace{1}=image;                                                 %base_image��Ϊ�߶ȿռ�ĵ�һ��ͼ��
%��ȡ�˲�������
h=[-1,0,1;-2,0,2;-1,0,1];      % soble ����˲�ģ��

%����ÿ��ĳ߶�
sigma=zeros(1,layers);
for i=1:1:layers
    sigma(i)=sigma_1*ratio^(i-1);%ÿ��ĳ߶�
end

%% ���������Գ߶ȿռ�
for i=2:1:layers
    %֮ǰ��ķ�������ɢ��ĵ�ͼ��,�����ݶ�֮ǰ����ƽ����Ŀ����Ϊ����������
    prev_image=Nonelinear_Scalespace{i-1};
    prev_image=imresize(prev_image,1/scale_value,'bilinear');
    windows_size=2*round(2*sigma_2)+1;
    W=fspecial('gaussian',[windows_size,windows_size],sigma_2);
    prev_smooth=imfilter(prev_image,W,'replicate');
    
    %����֮ǰ�㱻ƽ��ͼ���x��y�����һ���ݶ�
    Lx=imfilter(prev_smooth,h ,'replicate');
    Ly=imfilter(prev_smooth,h','replicate');   
    
    %ÿ�ε���ʱ����Ҫ���¶Աȶ�����k
    [k_percentile]=K_percentile_value(Lx,Ly,perc);
    if(which_diff==1)
        [diff_c]=pm_g1(Lx,Ly,k_percentile);
    elseif(which_diff==2)
        [diff_c]=pm_g2(Lx,Ly,k_percentile);
    else
        [diff_c]=weickert_diffusivity(Lx,Ly,k_percentile);
    end
    
    %���㵱ǰ��߶�ͼ��
    step=1/2*(sigma(i)^2-sigma(i-1)^2);%��������
    Nonelinear_Scalespace{i}=AOS(prev_image,step,diff_c);  %nonelinear_space: prev_image��ʾ֮ǰ��ķ�������ɢ��ĵ�ͼ�� step��ʾ������diff_c��ʾ��ɢϡ��C��
end
end

%% ��ɢϵ�����㺯��1
function [g1]=pm_g1(Lx,Ly,k)
%�ú�������PM����ϵ��g1,Lx��ˮƽ����ĵ�����Ly����ֱ����ĵ���
%k��һ���Աȶ����Ӳ�����k��ȡֵһ�����ͳ������
%g1=exp(-(Lx^2+Ly^2)/k^2)

g1=exp(-(Lx.^2+Ly.^2)/k^2);

end

%% ��ɢϵ�����㺯��2
function [g2]=pm_g2(Lx,Ly,k)
%�ú�������PM���̵���ɢϵ�����ڶ��ַ���
%Lx��Ly�ֱ���ˮƽ�������ֱ����Ĳ�֣�k�ǶԱȶ����Ӳ���
%g2=1/(1+(Lx^2+Ly^2)/(k^2)),����kֵ��ȷ��һ����ͨ��ͳ�Ʒ����õ�

g2=1./(1+(Lx.^2+Ly.^2)/(k^2));
end

%% ��ɢϵ�����㺯��3
function [g3]=weickert_diffusivity(Lx,Ly,k)
%�����������weickert����ϵ��
%Lx��Ly��ˮƽ�������ֱ�����һ�ײ���ݶȣ�k�ǶԱȶ�ϵ��
%k��ȡֵһ��ͨ��ͳ�Ʒ����õ�
%g3=1-exp(-3.315/((Lx^2+Ly^2)/k^4))

g3=1-exp(-3.315./((Lx.^2+Ly.^2).^4/k^8));
end

%%  Kֵ����
function [k_percentile]=K_percentile_value(gradient_x,gradient_y,perc)
%�ú�������һ���ԱȶȲ���k,����ԱȶȲ������ڼ�����ɢϵ��
%gradient_x��ˮƽ������ݶȣ�gradient_y����ֱ������ݶ�
%perc���ݶ�ֱ��ͼ�İٷ�λ����Ĭ��ȡֵ��0.7��k��ȡֵ��������ٷ�λ��ȷ��
%����ϵ����������k������������˶�����ͬ���ݶ�ֵ�����kֵ�ϴ��򴫵�ϵ��ֵ�ϴ�
%�����ɢ��ƽ�����أ���˿��Կ��������Ҫ����ϸ����Ҫ��С��kֵ
%�ú����Զ�����kֵ��������ﲻ��Ҫָ��bin�Ĵ�С

%ֱ��ͼ���
unit=0.005;

%���Ա߽�����ݶȵ����ֵ
gradient=sqrt(gradient_x.^2+gradient_y.^2);
[M,N]=size(gradient);
temp_gradient=gradient(2:M-1,2:N-1);

%���Ա߽����ֱ��ͼ
temp_gradient=temp_gradient(temp_gradient>0);
max_gradient=max(max(temp_gradient));
min_gradient=min(min(temp_gradient));
temp_gradient=round((temp_gradient-min_gradient)/unit);
nbin=round((max_gradient-min_gradient)/unit);
hist=zeros(1,nbin+1);
[M1,N1]=size(temp_gradient);
sum_pix=M1*N1;%���������ݶȸ���

%����ֱ��ͼ
for i=1:1:M1
    for j=1:1:N1
        hist(temp_gradient(i,j)+1)=hist(temp_gradient(i,j)+1)+1;
    end
end

%ֱ��ͼ�ٷ�λ
nthreshold=perc*sum_pix;
nelements=0;
temp_i=0;
for i=1:1:nbin+1
    nelements=nelements+hist(i);
    if(nelements>=nthreshold)
        temp_i=i;
        break;
    end
end
k_percentile=(temp_i-1)*unit+min_gradient;
end
