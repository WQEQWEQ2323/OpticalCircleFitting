clc;clear;close all;
% here, formula in filename means 1st hankel transform of box function

%% select input image
tic
[filename, pathname] = uigetfile({'*.png'; '*.jpg'; '*.gif'}, '选择图片');
% 没有图像
if isequal(filename, 0)
    return;
end
bg = imread(filename);
I=im2double(bg);
%I=rot90(I,1);
%I(I==0)=-1; phase image
%% parameter
f=0.1;
lambda=4.5e-7; %wavelength, meters
k=2*pi/lambda; %wavenumber k, m^-1
delta0=2e-6;c0=2000;L0=delta0*c0;%L0=4mm

%% 填充0 / enlarge FOV by Filling 0
I=F_ZerosFilling(I,c0,6*c0);
L=L0*6;c0=c0*6;delta=delta0; %L=2.4cm
x=(-c0/2:c0/2-1)*delta;[X,Y]=meshgrid(x,x);[theta,r]=cart2pol(X,Y);

%% 第一次传播 angular spectrum
df=lambda/L;f0=(-c0/2:c0/2-1)*df;
[AX,AY]=meshgrid(f0,f0);
C0=jfft2(I);
D0=C0.*exp(1i*k*f*sqrt(1-AX.^2-AY.^2)).*circ(sqrt(AX.^2+AY.^2));
%D0(isnan(D0))=0;
E1=jifft2(D0);
clear C0;clear D0;
%% 第二次传播
E1=E1.*exp(-1i*k/2/f*r.^2).*circ(r/(12e-3)); %radius of lens1 :1.2cm
C1=jfft2(E1);
D1=C1.*exp(1i*k*f*sqrt(1-AX.^2-AY.^2)).*circ(sqrt(AX.^2+AY.^2));
%D1(isnan(D1))=0;
E2=jifft2(D1);
clear C1;clear D1;
imagesc(x,x,abs(E2));colormap jet;colorbar;


%% u here u is gradient operator
R1=k/f*2e-6;R2=2*R1;%% R2不一定等于2倍R1 (the ratio of R2/R1 2 can change if the presentation can be improved)
x2=(-500:499)*delta;[X2,Y2]=meshgrid(x2,x2);[theta2,r2]=cart2pol(X2,Y2); % aperature of meta is 2mm(diameter) 2mm/delta = 1000pixel
u=(transmission_function(r2,R2)-transmission_function(r2,R1)).*exp(1i*theta2).*circ(r2/1e-3);
u=F_ZerosFilling(u,1000,c0);
u=u/max(max(u));
%% v
v=zeros(1000,1000);
for n=1:19 
    disp(strcat('n=',num2str(n)));
    R3=6*n*4e-6*k/f;%这里制造周期性+1,-1,+1,-1 (periodic reverse(pi phase) to introduce random/chaos for point offset from circle center)
    v=v+2*(-1)^n*(transmission_function(r2,R3));
end
n=20
R3=6*n*4e-6*k/f;
v=v+(-1)^n*(transmission_function(r2,R3));
v=v.*exp(-1i*theta2).*circ(r2/1e-3);

v=F_ZerosFilling(v,1000,c0);
v=v/max(max(v));

%% element
t=real(u.*v);
t=t/max(max(abs(t)));
%% double ring(simplified element)
kernel_1=ones(c0,c0);
kernel_1(r>0.994e-3)=0;
kernel_1(r<0.918e-3)=-1;
kernel_1(r<0.852e-3)=0;
%% 第三次传播
E2_new=E2.*t.*circ(r/1e-3); % aperature: 2mm
C2=jfft2(E2_new);
D2=C2.*exp(1i*k*f*sqrt(1-AX.^2-AY.^2)).*circ(sqrt(AX.^2+AY.^2));
E3=jifft2(D2);
clear C2;clear D2;
%% 第四次传播
E3=E3.*exp(-1i*k/2/f*r.^2).*circ(r/12e-3);
C3=jfft2(E3);
D3=C3.*exp(1i*k*f*sqrt(1-AX.^2-AY.^2)).*circ(sqrt(AX.^2+AY.^2));
E4=jifft2(D3);
clear C3;clear D3; 

E5=F_UnZerosFilling(E4,2400); % why is 2400 instead of 2000? because we need same background for nine output images(nine sub-boxes) to compute efficiency
E5=abs(E5);
E5=E5/max(max(E5));
imagesc(x,x,E5);colormap jet;colorbar;
%imagesc(x,x,abs(E4));colormap jet;colorbar;

load chirp;
sound(y,Fs);
%%
T=abs(E4).^2;
T=F_UnZerosFilling(T,2400);
T=255*T/max(T(:));
%T=T/max(max(T));
bg = F_UnZerosFilling(flip(flip(I,2),1),2400);
nine_square = cell(9,1);
nine_square{1}=[1,1,800,800];
nine_square{2}=[801,1,1600,800];
nine_square{3}=[1601,1,2400,800];
nine_square{4}=[1,801,800,1600];
nine_square{5}=[801,801,1600,1600];
nine_square{6}=[1601,801,2400,1600];
nine_square{7}=[1,1601,800,2400];
nine_square{8}=[801,1601,1600,2400];
nine_square{9}=[1601,1601,2400,2400];
focusing_efficiency = zeros(9,1);
figure;
imagesc(T);
hold on;
for i = 1:9
        poi = nine_square{i};
        x1 = poi(1);y1 = poi(2);x2 = poi(3);y2 = poi(4);
        sub_bg_mask = bg(y1:y2, x1:x2);
        subImg = T(y1:y2, x1:x2);
        %[best_loGcal_sum, peak_mask] = generateLocalSum(subImg);
        row = y2-y1+1;col = x2-x1+1;
        [Y,X]=meshgrid(1:col,1:row);
        Radial_length = sqrt((X-(1+row)/2).^2+(Y-(1+col)/2).^2);
        peak_mask_A = Radial_length  < 19/2;
    
        [best_local_sum, peak_mask] = generateLocalSum(subImg,sub_bg_mask);
    
        % 计算方框内像素的最大灰度值
        efficiency(i) = computeLocalAverageAroundPeak(subImg,peak_mask);
    
        % 在图像上绘制方框并显示最大灰度值
        rectangle('Position', [x1,y1,x2-x1,y2-y1], 'EdgeColor', 'r');
        text((x1+x2)/2-40, y1 + 40, sprintf('%.5f', efficiency(i)), 'Color', 'r', 'FontSize', 20, 'FontWeight', 'bold'); % 使用比值，无量纲，更贴近cv中的识别问题,scale/size independent
        text((x1+x2)/2-40, y1 + 120, sprintf('%.2f', max(subImg(:))), 'Color', 'g', 'FontSize', 20, 'FontWeight', 'bold'); % 我实际上更喜欢这个，能够说明问题，找到核质比更大的细胞

        hold on;
    end
%% old metric
%{
AA=T(1:400,1:400);aa=max(max(AA));
BB=T(1:400,801:1200);bb=max(max(BB));
CC=T(1:400,1601:2400);cc=max(max(CC));
DD=T(801:1200,1:400);dd=max(max(DD));
EE=T(801:1200,801:1200);ee=max(max(EE));
FF=T(801:1200,1601:2400);ff=max(max(FF));
GG=T(1601:2400,1:400);gg=max(max(GG));
HH=T(1601:2400,801:1200);hh=max(max(HH));
II=T(1601:2400,1601:2400);ii=max(max(II));
xx=round([aa,bb,cc,dd,ee,ff,gg,hh,ii]*255)
yy=sort(xx,'descend');
if aa==yy(1)
    yy(1)/yy(2)
else
    aa/yy(1)
end
toc
imagesc(T);colormap jet ;axis equal;hold on;
%}

function [best_local_sum, best_mask] = generateLocalSum(input_img, input_bg_mask, radius)
    if nargin < 3
        radius = 19/2;
    end

    if size(input_img,3) == 3
        input_img = rgb2gray(input_img);
    end
    input_img(~input_bg_mask) = 0; % actual useless for simulation because of noiseless
    [rows, cols] = size(input_img);

    % Step 1: 找所有最大值位置
    max_val = max(input_img(:));
    [max_y_list, max_x_list] = find(input_img == max_val);

    % 初始化
    best_local_sum = 0;
    best_mask = false(rows, cols);

    % 遍历所有最大值点
    for k = 1:length(max_y_list)
        y0 = max_y_list(k);
        x0 = max_x_list(k);

        % Step 2: 构建 mask
        [X, Y] = meshgrid(1:cols, 1:rows);
        R = sqrt((X - x0).^2 + (Y - y0).^2);
        mask = R <= radius;

        % Step 3: 计算局部和
        local_sum = sum(input_img(mask));

        % Step 4: 找最大 local_sum 的 mask
        if local_sum > best_local_sum
            best_local_sum = local_sum;
            best_mask = mask;
        end
    end
end


%% 最高峰周围的能量比上总能量
function focusing_efficiency = computeLocalAverageAroundPeak(input_img, peak_mask)
    % Step 1: 计算区域内灰度总和与输入图片总灰度
    local_sum = sum(input_img(peak_mask));
    global_sum = sum(input_img(:));

    % Step 2: 灰度之比
    focusing_efficiency = local_sum / global_sum;
    fprintf('主峰周围总灰度 = %.2f 虚线区域总灰度 = %.2f，比值 = %.4f\n', ...
        local_sum, global_sum, focusing_efficiency);
end