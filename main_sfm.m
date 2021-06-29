%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this toolbox is an addition to the toolbox provided by the authors of
% MLPnP, CEPPnP and OPnP
% we extended it to show the use of EKFPnP
%
% Copyright (C) <2018>  <MA.Mehralian>
%
%     email: ma.mehralian@gmail.com
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_sfm()
clear; clc;
IniToolbox;

warning off;
if ~exist('vl_sift', 'file')
    if ~exist('vl_sift', 'file')
       disp('"vlfeat/toolbox/vl_setup.m" not found. Please download "vlfeat" and put it next to this code"')
    end
    run('vlfeat/toolbox/vl_setup');
end

ekf_err=[]; ceppnp_err=[]; mlpnp_err=[];
dataset = 'dino'; %'temple'; %
for ex=1:100
    [data, K] = read_data_middlebury(dataset);
    invK = inv(K);
    K_f = K(1,1);
    K_c = K(1:2,3);
    
    %----- initial reconstruction using groundtruth
    n = length(data{1}.matches);
    x1 = data{1}.f(:, data{1}.matches(1,:));
    x2 = data{2}.f(:, data{1}.matches(2,:));
    x1_= invK * x1;
    x2_= invK * x2;
    P = cat(3, data{1}.proj, data{2}.proj);
    X = triangulate(P, cat(3, x1_, x2_));
    % figure; structure_plot(P, X);
    X_f = nan(3, length(data{1}.f));
    X_f(:,data{1}.matches(1,:)) = X(1:3,:);
    
    
    %---- filter initialization
    P1_q = p2state( data{1}.proj );
    P2_q = p2state( data{2}.proj );
    f = ekfpnp(P1_q, P2_q, [K(1,1) K(1:2,3)']);
    
    ekf_P = P;
    ceppnp_P = P;
    mlpnp_P = P;
    grt_P = P;
    
    data{1}.ekf_X = X_f;
    data{1}.ceppnp_X = X_f;
    data{1}.mlpnp_X = X_f;
    
    
    %---- PROCESS SATRT
    s_len = length(data)-1;
    for step=2:s_len
        %--- find common feature with last 3 images which have 3D correspond
        [ids, ia, ib] = intersect(data{step-1}.matches(2,:), data{step}.matches(1,:));
        fprintf('STEP %g: %g measurments\n', step, length(ia));
        
        f0_idx = data{step-1}.matches(1,ia);
        f1_idx = data{step}.matches(1,ib); %%OR%% data{step-1}.matches(2,ia)  %%OR%% ids
        f2_idx = data{step}.matches(2,ib);
        
        x0 = data{step-1}.f(:, f0_idx);
        x1 = data{step}.f(:, f1_idx);
        x2 = data{step+1}.f(:, f2_idx);
        %--- EKFPnP
        X_ekf = data{step-1}.ekf_X(1:3, f0_idx);
        z = reshape(x2(1:2,:),[],1);
        f = f.step_simple(X_ekf, z);
        ekf_P(:,:,end+1) = state2p( f.s_k_k(1:7) );
        
        %----------
        npt = length(ia);
        x_c = (x2(1:2,:)-repmat(K_c,1,npt))/K_f;
        homx = [x_c; ones(1,npt)];
        v = normc(homx);
        Cu = zeros(2,2,npt);
        Evv = zeros(3,3,npt);
        cov = zeros(9,size(Cu,3));
        nl = 1;
        for id = 1:npt
            Cu(:,:,id) = diag([nl^2 nl^2]);
            cov_proj = K\diag([nl^2 nl^2 0])/K';
            J = (eye(3)-(v(:,id)*v(:,id)')/(v(:,id)'*v(:,id)))/norm(homx(:,id));
            Evv(:,:,id) = J*cov_proj*J';
            cov(:,id) = reshape(Evv(:,:,id),9,1);
        end
        %--- CEPPNP
        X_ceppnp = data{step-1}.ceppnp_X(1:3, f0_idx);
        mX = X_ceppnp - repmat(mean(X_ceppnp,2),1,size(X_ceppnp,2));
        [R_ceppnp,t_ceppnp]=CEPPnP(mX, x_c, Cu);
        t_ceppnp = t_ceppnp - R_ceppnp * mean(X_ceppnp,2);
        ceppnp_P(:,:,end+1) = [R_ceppnp, t_ceppnp];
        
        %--- MLPnP
        X_mlpnp = data{step-1}.mlpnp_X(1:3, f0_idx);
        [R_mlpnp,t_mlpnp]=MLPNP_without_COV(X_mlpnp, v, []);
        mlpnp_P(:,:,end+1) = [R_mlpnp, t_mlpnp];
        
        grt_P(:,:,end+1) = data{step+1}.proj;
        
        
        %--- find 3D corespond of new points
        m = size(data{step}.matches, 2);
        idx = true(1,m);
        idx(ib) = 0;
        ibn = find(idx);
        f1_idx_new = data{step}.matches(1,ibn);
        f2_idx_new = data{step}.matches(2,ibn);
        x1_= invK * data{step}.f(:, f1_idx_new);
        x2_= invK * data{step+1}.f(:, f2_idx_new);
        
        %--- ekf triangulation
        P = ekf_P(:,:,end-1:end);
        X_ekf_new = triangulate(P, cat(3, x1_, x2_));
        data{step}.ekf_X = nan(3, length(data{step}.f));
        data{step}.ekf_X(:, f1_idx) = data{step-1}.ekf_X(1:3, f0_idx);
        data{step}.ekf_X(:, f1_idx_new) = X_ekf_new(1:3,:);
        
        %--- ceppnp triangulation
        P = ceppnp_P(:,:,end-1:end);
        X_ceppnp_new = triangulate(P, cat(3, x1_, x2_));
        data{step}.ceppnp_X = nan(3, length(data{step}.f));
        data{step}.ceppnp_X(:, f1_idx) = data{step-1}.ceppnp_X(1:3, f0_idx);
        data{step}.ceppnp_X(:, f1_idx_new) = X_ceppnp_new(1:3,:);
        
        %--- mlpnp triangulation
        P = mlpnp_P(:,:,end-1:end);
        X_mlpnp_new = triangulate(P, cat(3, x1_, x2_));
        data{step}.mlpnp_X = nan(3, length(data{step}.f));
        data{step}.mlpnp_X(:, f1_idx) = data{step-1}.mlpnp_X(1:3, f0_idx);
        data{step}.mlpnp_X(:, f1_idx_new) = X_mlpnp_new(1:3,:);
    end
    [ekf_err(:,:,ex), ceppnp_err(:,:,ex), mlpnp_err(:,:,ex)] = getErrors(p2state(grt_P), p2state(ekf_P), p2state(ceppnp_P), p2state(mlpnp_P));
end

%--- plot error results
ekf_err = mean(ekf_err,3);
ceppnp_err = mean(ceppnp_err,3);
mlpnp_err = mean(mlpnp_err,3);
plotErros(ekf_err, ceppnp_err, mlpnp_err);
saveas(gcf, ['real_expriment_sfm_',dataset,'.eps'], 'epsc')

%--- plot SFM results
X_ekf = [];  X_ceppnp = [];  X_mlpnp = [];
for i=1:s_len
    ix_ekf = ~isnan(data{i}.ekf_X(1,:));
    ix_ceppnp = ~isnan(data{i}.ceppnp_X(1,:));
    ix_mlpnp = ~isnan(data{i}.mlpnp_X(1,:));
    X_ekf = [X_ekf data{i}.ekf_X(:,ix_ekf)];
    X_ceppnp = [X_ceppnp data{i}.ceppnp_X(:,ix_ceppnp)];
    X_mlpnp = [X_mlpnp data{i}.mlpnp_X(:,ix_mlpnp)];
end

figure; hold on; axis equal;
structure_plot(grt_P, nan(3,1), 'cam_color', 'k', 'cam_fill', 'g');
structure_plot(ekf_P, X_ekf, 'cam_color', 'r');
view(-180, -80);
title('EKFPnP');
saveas(gcf, ['real_expriment_sfm_',dataset,'_EKFPnP.eps'],'epsc')

% figure; hold on; axis equal;
% structure_plot(grt_P, nan(3,1), 'cam_color', 'k', 'cam_fill', 'g');
% structure_plot(ceppnp_P, X_ceppnp, 'cam_color', 'r');
% view(-180, -80);
% title('CEPPnP');
% saveas(gcf, ['real_expriment_sfm_',dataset,'_CEPPnP.eps'],'epsc')

figure; hold on; axis equal;
structure_plot(grt_P, nan(3,1), 'cam_color', 'k', 'cam_fill', 'g');
structure_plot(mlpnp_P, X_mlpnp, 'cam_color', 'r');
view(-180, -80);
title('MLPnP');
saveas(gcf, ['real_expriment_sfm_',dataset,'_MLPnP.eps'],'epsc')
end
%%
function name = get_file_name(path)
tmp = strsplit(path,'/');
name = tmp{end};
end
%%
function data = match_features(images, SCALE, KNN_RATION)
% f = fetures
% d = decriptors

PEAK_THRESH = 0.001;
LEVELS = 5;
data = [];

% ----- feature extraction
for i=1:length(images)
    img = im2single(rgb2gray(imread(images{i})));
    img = imresize(img, SCALE);
    [data{i}.f, data{i}.d] = vl_sift(img);%,'PeakThresh', PEAK_THRESH, 'Levels', LEVELS);
    disp([num2str(length(data{i}.f)), ' features are extracted from ', get_file_name(images{i})]);
    %{
    figure; imshow(img);
    hold on;
    scatter(data{i}.f(1,:), data{i}.f(2,:),100,'.r');
    %}
end
    
% ----- match features using 2nd KNN method and RANSAC
disp('Feature matching: ');
for i = 1:length(images)-1
    j = i+1;
    fprintf('Feature matching (%s to %s) --------\n', get_file_name(images{i}), get_file_name(images{j}));
    % --- 2nd KNN + cross check
    [knn_ix_ij,knn_d_ij] = knnsearch(data{i}.d', data{j}.d', 'k',2);
    knn_ratio_ij = knn_d_ij(:,1)./knn_d_ij(:,2);
    bst_ix_ij = find(knn_ratio_ij<KNN_RATION);
    
    [knn_ix_ji,knn_d_ji] = knnsearch(data{j}.d', data{i}.d', 'k',2);
    knn_ratio_ji = knn_d_ji(:,1)./knn_d_ji(:,2);
    bst_ix_ji = find(knn_ratio_ji<KNN_RATION);
    
    m_1 = sortrows([knn_ix_ij(bst_ix_ij,1) bst_ix_ij]);
    m_2 = sortrows([bst_ix_ji knn_ix_ji(bst_ix_ji,1)]);
     
    %----
    data{i}.matches = intersect(m_1, m_2, 'row')';
    fprintf('\t2nd KNN reduces features from %d to %d\n', length(data{i}.d), ...
        length(data{i}.matches));
    
    % --- RANSAC
    p_1 = data{i}.f(1:2,data{i}.matches(1,:));
    p_2 = data{j}.f(1:2,data{i}.matches(2,:));
    p_1(3,:) = 1;     p_2(3,:)= 1;
    [~,inliers] = fundamental_projective(p_1, p_2, 'use_ransac'); %FRansac(p_1, p_2);
    fprintf('\tRANSAC reduces features from %d to %d\n', ...
        length(data{i}.matches), sum(inliers));
    data{i}.matches = data{i}.matches(:,inliers);
    
    %--- plot
    %{
        img_i = imresize(im2single(rgb2gray(imread(images{i}))), SCALE);
        img_j = imresize(im2single(rgb2gray(imread(images{j}))), SCALE);
        f1 = data{i}.f(1:2, data{i}.matches(1,1:10:end));
        f2 = data{j}.f(1:2, data{i}.matches(2,1:10:end));
        match_plot(img_i, img_j, f1, f2);
    %}
end
    
end
%%
function [data, K] = read_data_ycb()
% http://www.ycbbenchmarks.com/
% http://ycb-benchmarks.s3-website-us-east-1.amazonaws.com/
K = hdf5read('data/035_power_drill_berkeley_rgb_highres/035_power_drill/calibration.h5','/N1_rgb_K')';
data_dir = 'data/035_power_drill_berkeley_rgb_highres/035_power_drill/*.jpg';
pose_dir = 'data/035_power_drill_berkeley_rgb_highres/035_power_drill/poses/*.h5';


images = dir(data_dir);
images = arrayfun(@(f)[f.folder,'/', f.name], images, 'UniformOutput', false);
poses = dir(pose_dir);
poses = arrayfun(@(f)[f.folder,'/', f.name], poses, 'UniformOutput', false);

end

%%
function [data, K] = read_data_epfl(DATASET)
% dataset from https://icwww.epfl.ch/multiview/denseMVS.html

SCALE = 0.25;
KNN_RATION = 0.7;
PEAK_THRESH = 0.001;
LEVELS = 5;

K = [ SCALE*2759.48 0 SCALE*1520.69 
0 SCALE*2764.16 SCALE*1006.81 
0 0 1 ];
% if exist(['data_',DATASET,'_',num2str(SCALE),'_',num2str(KNN_RATION),'.mat'])
%     load(['data_',DATASET,'_',num2str(SCALE),'_',num2str(KNN_RATION),'.mat']);
% else
    images = dir(['data/', DATASET, '/*.png']);
    images = arrayfun(@(f)[f.folder,'/', f.name], images, 'UniformOutput', false);
    cams = dir(['data/', DATASET, '/*.camera']);
    cams = arrayfun(@(f)[f.folder,'/', f.name], cams, 'UniformOutput', false);

    data = match_features(images, SCALE, KNN_RATION);
    tmp = what(['data/', DATASET, '/']);
    for i=1:length(images)
        data{i}.img_name = strrep(images{i}, [tmp.path,'/'],''); 
        %--- homogen features
        f = [data{i}.f(1:2,:); ones(1,size(data{i}.f,2))];
        %f = inv(K) * f;
        %f = f ./ f(3,:);
        data{i}.f = f;
        
        %--- cam pose
        cam = importdata(cams{i});
        data{i}.pose = [cam(5:7, :) cam(8, :)'];
        data{i}.proj = pInv(data{i}.pose);
    end
    %save(['data_',DATASET,'_',num2str(SCALE),'_',num2str(KNN_RATION),'.mat'], 'data');
% end

%{
Ps = cellfun(@(x) x.proj, data, 'UniformOutput', false);
Ps = cat(3, Ps{:});
figure; structure_plot(Ps, [0,0,0]');
%}
end
%%
function [data, K] = read_data_middlebury(DATASET)
% dataset from http://vision.middlebury.edu/mview/data/
if (strcmp(DATASET, "temple"))
    data_dir = 'data/templeRing/*.png';
    cam_file = 'data/templeRing/templeR_par.txt';
    cam_range = 13:31;
elseif(strcmp(DATASET, "dino"))
    data_dir = 'data/dino/*.png';
    cam_file = 'data/dino/dino_par.txt';
    cam_range = 1:48;
        
else
	disp('Unknown middlebury dataset')
end


KNN_RATION = 0.6;

K =[1520.4  0      302.3
    0       1525.9 246.9
    0       0      1];

images = dir(data_dir);
images = arrayfun(@(f)[f.folder,'/', f.name], images, 'UniformOutput', false);
cams = importdata(cam_file, ' ', 1);

idx = cam_range;
%images = images(idx);
cams.data = cams.data(idx,:);
cams.textdata = cams.textdata(idx+1,:);
data = match_features(images, 1, KNN_RATION);
for i=1:length(images)
    %--- homogen features
    f = [data{i}.f(1:2,:); ones(1,size(data{i}.f,2))];
    %f = inv(K) * f;
    %f = f ./ f(3,:);
    data{i}.f = f;
    
    %--- cam pose      cams.textdata(i+1)
    K = reshape(cams.data(i,1:9),[3,3])';
    R = reshape(cams.data(i,10:18),[3,3])';
    t = cams.data(i,19:21)';
    data{i}.proj = [R t];
    data{i}.pose = pInv(data{i}.proj);
end
%save(['data_',DATASET,'_',num2str(KNN_RATION),'.mat'], 'data');

%{
Ps = cellfun(@(x) x.proj, data, 'UniformOutput', false);
Ps = cat(3, Ps{:});
figure; structure_plot(Ps, nan(3,1));
%}
end
%%
function [ekf_err, ceppnp_err, mlpnp_err]=getErrors(grt_x, ekf_x, ceppnp_x, mlpnp_x)
%{
a=sqrt(sum((ekf_x(1:3,:)-grt_x(1:3,:)).^2))
b=sqrt(sum((pnp_x(1:3,:)-grt_x(1:3,:)).^2))
figure; plot(a); hold on; plot(b);
%}

for i=1:size(grt_x,2)
    %--- transition errors
    mt_ekf(i) = norm(ekf_x(1:3,i)-grt_x(1:3,i))./norm(grt_x(1:3,i))*100;
    mt_ceppnp(i) = norm(ceppnp_x(1:3,i)-grt_x(1:3,i))./norm(grt_x(1:3,i))*100;
    mt_mlpnp(i) = norm(mlpnp_x(1:3,i)-grt_x(1:3,i))./norm(grt_x(1:3,i))*100;
    
    %--- rotatin error
    ekf_dq = quatmultiply(quatinv(grt_x(4:7,i)'), ekf_x(4:7,i)');
    [ekf_dq_yaw, ekf_dq_pitch, ekf_dq_roll] = quat2angle(ekf_dq);
    mr_ekf(i) = rad2deg(norm([ekf_dq_yaw, ekf_dq_pitch, ekf_dq_roll]));
    
    ceppnp_dq = quatmultiply(quatinv(grt_x(4:7,i)'), ceppnp_x(4:7,i)');
    [ceppnp_dq_yaw, ceppnp_dq_pitch, ceppnp_dq_roll] = quat2angle(ceppnp_dq);
    mr_ceppnp(i) = rad2deg(norm([ceppnp_dq_yaw, ceppnp_dq_pitch, ceppnp_dq_roll]));

    mlpnp_dq = quatmultiply(quatinv(grt_x(4:7,i)'), mlpnp_x(4:7,i)');
    [mlpnp_dq_yaw, mlpnp_dq_pitch, mlpnp_dq_roll] = quat2angle(mlpnp_dq);
    mr_mlpnp(i) = rad2deg(norm([mlpnp_dq_yaw, mlpnp_dq_pitch, mlpnp_dq_roll]));

end
ekf_err=[mt_ekf; mr_ekf];
ceppnp_err=[mt_ceppnp; mr_ceppnp];
mlpnp_err=[mt_mlpnp; mr_mlpnp];
end
%%
function plotErros(ekf_err, ceppnp_err, mlpnp_err)
figure;
subplot(2,1,1);
hold on;
grid on;
plot(ekf_err(1,:), '-r', 'linewidth', 1);
plot(ceppnp_err(1,:), '-b', 'linewidth', 1);
plot(mlpnp_err(1,:), '-k', 'linewidth', 1);
title('Mean Translation Error');
xlabel('n');
ylabel('Translation Error (%)');
legend({'EKFPnP', 'MLPnP'}, 'Location','northwest');
legend('EKFPnP', 'CEPnP', 'MLPnP');
ylim([0, 100])

subplot(2,1,2);
hold on;
grid on;
plot(ekf_err(2,:), '-r', 'linewidth', 1);
plot(ceppnp_err(2,:), '-b', 'linewidth', 1);
plot(mlpnp_err(2,:), '-k', 'linewidth', 1);
title('Mean Rotation Error');
xlabel('n');
ylabel('Rotation Error (degrees)');
legend({'EKFPnP', 'MLPnP'}, 'Location','northwest');
legend('EKFPnP', 'CEPnP', 'MLPnP');
end