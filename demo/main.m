clc;clear;close all;


addpath('.\src');
addpath('.\extra');
load('.\data\grace4.mat');

gt=1; gt_threshold = 3;
if size(I1,3)>1
img1 = I1; img2 = I2; GT = H;
else
    img1(:,:,1) = I1;img1(:,:,2) = I1; img1(:,:,3) = I1;
    img2(:,:,1) = I2;img2(:,:,2) = I2; img2(:,:,3) = I2;
    GT = H;
end
[tar_feat,tar_desc, ref_feat, ref_desc] = sift_process(img1,img2);

newfig=zeros(max(size(img1,1),size(img2,1)), size(img1,2)+size(img2,2),3);
newfig(:,1:size(img1,2),:) = img1;
newfig(1:size(img2,1),(size(img1,2)+1):end,:)=img2;
newfig=uint8(newfig);
figure;
image(newfig);title('original image pair');
axis image;axis off;

figure;
image(newfig);title('feature extraction');
axis image; axis off;
f1 = tar_feat; f2 = ref_feat;
f2(1,:) = ref_feat(1,:)+size(img1,2);
h1 = vl_plotframe(f1);
h2 = vl_plotframe(f2);
set(h1,'color','g','linewidth',2);
set(h2,'color','r','linewidth',2);
hold off

t_xd = 30; t_yd = 80;

%% Progressive smoothness consensus 
tic;
f_matches = PSC(tar_feat, ref_feat, tar_desc, ref_desc);
toc;

figure;
image(newfig);title('PSC');
axis image; axis off;
f1 = tar_feat; f2 = ref_feat;
f2(1,:) = ref_feat(1,:)+size(img1,2);hold on;
ind = 1:size(f_matches,2);
if length(ind)>200
    pind = randperm(numel(ind),200);
else
    pind = ind;
end

plot_matches = f_matches(:,pind);
plot(f1(1,plot_matches(1,:)), f1(2,plot_matches(1,:)), 'y.', f2(1,plot_matches(2,:)), f2(2,plot_matches(2,:)), 'y.', 'MarkerSize', 10.0);
hold on

if gt
    inliers = ground_truth_verification( tar_feat, ref_feat, plot_matches, GT, gt_threshold );
    line([f1(1,plot_matches(1,inliers));f2(1,plot_matches(2,inliers))],[f1(2,plot_matches(1,inliers));f2(2,plot_matches(2,inliers))],'linewidth',1,'color','b');
    hold on;
    line([f1(1,plot_matches(1,~inliers));f2(1,plot_matches(2,~inliers))],[f1(2,plot_matches(1,~inliers));f2(2,plot_matches(2,~inliers))],'linewidth',1,'color','r');
    
    inliers = ground_truth_verification( tar_feat, ref_feat, f_matches, GT, gt_threshold );
    psc= sprintf('#matches: %d, #inliers: %d, PC: %.4f', size(f_matches,2), sum(inliers), sum(inliers)/size(f_matches, 2))
    text(t_xd, t_yd, psc, 'FontUnits', 'pixels', 'FontSize', 10, 'Color', [0.95,0.95,0.95], 'BackgroundColor', [0.2,0.2,0.2]);  
else
    line([f1(1,f_matches(1,pind));f2(1,f_matches(2,pind))],[f1(2,f_matches(1,pind));f2(2,f_matches(2,pind))],'linewidth',1,'color','m');
    hold off;
end

%% Two-step strategy - CRC
tic;
sift_theta = 1.5;
[matches, ~] = vl_ubcmatch(tar_desc, ref_desc, sift_theta);
X = tar_feat(1:2,matches(1,:))';
Y = ref_feat(1:2,matches(2,:))';
% mismatch removal

[Xn, Yn] = normr(X, Y);
conf4 = CRC_init([]);
[indx, ~, ~] = CRC(Xn, Yn, conf4);

toc;
% plot result
figure;
image(newfig);title('CRC');
axis image; axis off;
f1 = tar_feat; f2 = ref_feat;
f2(1,:) = ref_feat(1,:)+size(img1,2);hold on;
f_matches = matches(:,indx);
ind = 1:size(f_matches,2);
if length(ind)>200
    pind = randperm(numel(ind),200);
else
    pind = ind;
end
plot_matches = f_matches(:,pind);
plot(f1(1,plot_matches(1,:)), f1(2,plot_matches(1,:)), 'y.', f2(1,plot_matches(2,:)), f2(2,plot_matches(2,:)), 'y.', 'MarkerSize', 10.0);
hold on

if gt
    inliers = ground_truth_verification( tar_feat, ref_feat, plot_matches, GT, gt_threshold );
    line([f1(1,plot_matches(1,inliers));f2(1,plot_matches(2,inliers))],[f1(2,plot_matches(1,inliers));f2(2,plot_matches(2,inliers))],'linewidth',1,'color','b');
    hold on;
    line([f1(1,plot_matches(1,~inliers));f2(1,plot_matches(2,~inliers))],[f1(2,plot_matches(1,~inliers));f2(2,plot_matches(2,~inliers))],'linewidth',1,'color','r');
    
    inliers = ground_truth_verification( tar_feat, ref_feat, f_matches, GT, gt_threshold );
     % sprintf('#features: %d, #matches: %d, #inliers: %d, ratio: %.4f', sum(size(f1,2)+size(f2,2)), size(f_matches,2), sum(inliers), sum(inliers)/size(f_matches, 2))
    crc = sprintf('#matches: %d, #inliers: %d, PC: %.4f', size(f_matches,2), sum(inliers), sum(inliers)/size(f_matches, 2))
    text(t_xd, t_yd, crc, 'FontUnits', 'pixels', 'FontSize', 10, 'Color', [0.95,0.95,0.95], 'BackgroundColor', [0.2,0.2,0.2]);  
else
    line([f1(1,f_matches(1,pind));f2(1,f_matches(2,pind))],[f1(2,f_matches(1,pind));f2(2,f_matches(2,pind))],'linewidth',1,'color','m');
    hold off;
end
