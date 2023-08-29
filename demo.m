clc;
clear;
close all;
%% 1 Import and display reference and image to be registered
file_image= 'G:\1\';
[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select Image',file_image);
image1=imread(strcat(pathname,filename));
[filename,pathname]=uigetfile({'*.*','All Files(*.*)'},'Select Image',file_image);
image2=imread(strcat(pathname,filename));

%% 2  Setting of initial parameters 
% Key parameters:
K_weight=3;                        % 各向异性力矩图的加权值，处于（1~10），默认设置：3
Max=3;                               % Number of levels in scale space，尺度空间的层数，默认设置：3
threshold = 0.4;                  % 特征点提取阈值，对SAR影像/强度图配色时，设置为：0.3；一般默认设置为：0.4
scale_value=2;                  % 尺度缩放比例值，默认设置：1.6
Path_Block=42;                   % 描述子邻域窗口大小， 默认设置：42；当需要更多特征点时，可以调大窗口。
K =0.15;   
sigma_1=1.6;                       %The first level of scale, 默认值是：1.6
ratio=2^(1/3);                     % scale ratio;
Scale ='NO';
%% 3 各向异性尺度空间
t1=clock;
disp('please waiting...');
tic;
[nonelinear_space_1]=nonelinear_space(image1,Max,scale_value);
[nonelinear_space_2]=nonelinear_space(image2,Max,scale_value);
disp(['构造各向异性尺度空间花费时间：',num2str(toc),'秒']);

%% 4  构建加权各向异性力矩图和相位一致性梯度及方向 
tic;
[harris_function_1,gradient_1,angle_1]=Gradient_Feature(nonelinear_space_1,Max,K_weight);
[harris_function_2,gradient_2,angle_2]=Gradient_Feature(nonelinear_space_2,Max,K_weight);
disp(['构建各向异性加权力矩图:',num2str(toc),'S']);

%% 5  feature point extraction
tic;
position_1=Harris_extreme(harris_function_1,gradient_1,angle_1,Max,threshold);
position_2=Harris_extreme(harris_function_2,gradient_2,angle_2,Max,threshold);
figure(1),imshow(image1); hold on; plot(position_1(:,1),position_1(:,2),'*r'); pause(0.01)
figure(2),imshow(image2); hold on; plot(position_2(:,1),position_2(:,2),'+g'); pause(0.01)

disp(['特征点提取花费时间:  ',num2str(toc),' S']);

%% 6 Lop-Polar Descriptor Constrained by HAPCG
tic;
descriptors_1=GGLOH_descriptors(gradient_1,angle_1,position_1,Path_Block);                                     
descriptors_2=GGLOH_descriptors(gradient_2,angle_2,position_2,Path_Block);  
disp(['特征描述子花费时间:  ',num2str(toc),'S']); 

%% 7 Nearest matching    
disp('Nearest matching')
[indexPairs,~] = matchFeatures(descriptors_1.des,descriptors_2.des,'MaxRatio',1,'MatchThreshold', 10);
matchedPoints_1 = descriptors_1.locs(indexPairs(:, 1), :);
matchedPoints_2 = descriptors_2.locs(indexPairs(:, 2), :);
%% Outlier removal 
uni1=[matchedPoints_1(:,[1 2]),matchedPoints_2(:,[1 2])];
[~,i,~]=unique(uni1,'rows','first');
matchedPoints_1 = matchedPoints_1(sort(i)',:);
matchedPoints_2 = matchedPoints_2(sort(i)',:);

%% FSC
[H,rmse,cor1,cor2]=FSC(matchedPoints_1,matchedPoints_2,'affine',3);
[clearedPoints1,clearedPoints2] = BackProjection(cor1,cor2,scale_value);
P3=Finepoints(clearedPoints1(:,[1,2]),clearedPoints2(:,[1,2]));
P1_subpixel=clearedPoints1(:,[1,2]);
P2_subpixel=P3;
figure; showMatchedFeatures(image1, image2, P1_subpixel(:,[1,2]), P2_subpixel(:,[1,2]), 'montage');
disp('The correct number of corresponding points is： '); disp(size(P1_subpixel,1));
[T,RMSE]=FSC(P1_subpixel(:,[1,2]),P2_subpixel(:,[1,2]),'affine',3);
disp(['RMSE of Matching results: ',num2str(RMSE),'  像素']);
t2=clock;
disp(['算法匹配总共花费时间  :',num2str(etime(t2,t1)),' S']);     
disp('registration result')
image_fusion(image2,image1,double(T));





