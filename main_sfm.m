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

ekf_err=[]; cepnp_err=[]; mlpnp_err=[];
for ex=1:50
    [data, K] = read_data_Temple('templeRing');
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
    f = ekfp(P1_q, P2_q, [K(1,1) K(1:2,3)']);
    
    ekf_P = P;
    cepnp_P = P;
    mlpnp_P = P;
    grt_P = P;
    
    data{1}.ekf_X = X_f;
    data{1}.cepnp_X = X_f;
    data{1}.mlpnp_X = X_f;
    
    
    %---- PROCESS SATRT
    s_len = length(data)-1;
    for step=2:s_len
        %--- find common feature with last 3 images which have 3D correspond
        [ids, ia, ib] = intersect(data{step-1}.matches(2,:), data{step}.matches(1,:));
        fprintf('\nSTEP %g: %g measurments\n', step, length(ia));
        
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
        %--- CEPNP
        X_cepnp = data{step-1}.cepnp_X(1:3, f0_idx);
        mX = X_cepnp - repmat(mean(X_cepnp,2),1,size(X_cepnp,2));
        [R_cepnp,t_cepnp]=CEPPnP(mX, x_c, Cu);
        t_cepnp = t_cepnp - R_cepnp * mean(X_cepnp,2);
        cepnp_P(:,:,end+1) = [R_cepnp, t_cepnp];
        
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
        
        %--- cepnp triangulation
        P = cepnp_P(:,:,end-1:end);
        X_cepnp_new = triangulate(P, cat(3, x1_, x2_));
        data{step}.cepnp_X = nan(3, length(data{step}.f));
        data{step}.cepnp_X(:, f1_idx) = data{step-1}.cepnp_X(1:3, f0_idx);
        data{step}.cepnp_X(:, f1_idx_new) = X_cepnp_new(1:3,:);
        
        %--- mlpnp triangulation
        P = mlpnp_P(:,:,end-1:end);
        X_mlpnp_new = triangulate(P, cat(3, x1_, x2_));
        data{step}.mlpnp_X = nan(3, length(data{step}.f));
        data{step}.mlpnp_X(:, f1_idx) = data{step-1}.mlpnp_X(1:3, f0_idx);
        data{step}.mlpnp_X(:, f1_idx_new) = X_mlpnp_new(1:3,:);
    end
    [ekf_err(:,:,ex), cepnp_err(:,:,ex), mlpnp_err(:,:,ex)] = getErrors(p2state(grt_P), p2state(ekf_P), p2state(cepnp_P), p2state(mlpnp_P));
end

%--- plot error results
ekf_err = mean(ekf_err,3);
cepnp_err = mean(cepnp_err,3);
mlpnp_err = mean(mlpnp_err,3);
plotErros(ekf_err, cepnp_err, mlpnp_err);
saveas(gcf, 'real_expriment_sfm.eps', 'epsc')

%--- plot SFM results
X_ekf = [];  X_cepnp = [];  X_mlpnp = [];
for i=1:s_len
    ix_ekf = ~isnan(data{i}.ekf_X(1,:));
    ix_cepnp = ~isnan(data{i}.cepnp_X(1,:));
    ix_mlpnp = ~isnan(data{i}.mlpnp_X(1,:));
    X_ekf = [X_ekf data{i}.ekf_X(:,ix_ekf)];
    X_cepnp = [X_cepnp data{i}.cepnp_X(:,ix_cepnp)];
    X_mlpnp = [X_mlpnp data{i}.mlpnp_X(:,ix_mlpnp)];
end

figure; hold on; axis equal;
structure_plot(grt_P, nan(3,1), 'cam_color', 'k', 'cam_fill', 'g');
structure_plot(ekf_P, X_ekf, 'cam_color', 'r');
view(-180, -80);
title('EKFPnP');
saveas(gcf, 'real_expriment_sfm_EKFPnP.eps','epsc')

% figure; hold on; axis equal;
% structure_plot(grt_P, nan(3,1), 'cam_color', 'k', 'cam_fill', 'g');
% structure_plot(cepnp_P, X_cepnp, 'cam_color', 'r');
% view(-180, -80);
% title('CEPnP');
% saveas(gcf, 'real_expriment_sfm_CEPnP.eps')

figure; hold on; axis equal;
structure_plot(grt_P, nan(3,1), 'cam_color', 'k', 'cam_fill', 'g');
structure_plot(mlpnp_P, X_mlpnp, 'cam_color', 'r');
view(-180, -80);
title('MLPnP');
saveas(gcf, 'real_expriment_sfm_MLPnP.eps','epsc')
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
    disp(['Feature extraction (', get_file_name(images{i}), ')']);
    [data{i}.f, data{i}.d] = vl_sift(img);%,'PeakThresh', PEAK_THRESH, 'Levels', LEVELS);
    %{
    figure; imshow(img);
    hold on;
    plot(data{i}.f(1,:), data{i}.f(2,:),'.');
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
    fprintf('\t2nd KNN reduced features from %d to %d\n', length(data{i}.d), ...
        length(data{i}.matches));
    
    % --- RANSAC
    p_1 = data{i}.f(1:2,data{i}.matches(1,:));
    p_2 = data{j}.f(1:2,data{i}.matches(2,:));
    p_1(3,:) = 1;     p_2(3,:)= 1;
    [~,inliers] = fundamental_projective(p_1, p_2, 'use_ransac'); %FRansac(p_1, p_2);
    fprintf('\tRANSAC  reduced features from %d to %d\n', ...
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
function [data, K] = read_data_Temple(DATASET)
% dataset from http://vision.middlebury.edu/mview/data/

KNN_RATION = 0.6;

K =[1520.4  0      302.3
    0       1525.9 246.9
    0       0      1];

images = dir(['data/', DATASET, '/*.png']);
images = arrayfun(@(f)[f.folder,'/', f.name], images, 'UniformOutput', false);
cams = importdata(['data/', DATASET, '/templeR_par.txt'], ' ', 1);

idx = 13:31;
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
save(['data_',DATASET,'_',num2str(KNN_RATION),'.mat'], 'data');

%{
Ps = cellfun(@(x) x.proj, data, 'UniformOutput', false);
Ps = cat(3, Ps{:});
figure; structure_plot(Ps, nan(3,1));
%}
end
%%
function [ekf_err, cepnp_err, mlpnp_err]=getErrors(grt_x, ekf_x, cepnp_x, mlpnp_x)
%{
a=sqrt(sum((ekf_x(1:3,:)-grt_x(1:3,:)).^2))
b=sqrt(sum((pnp_x(1:3,:)-grt_x(1:3,:)).^2))
figure; plot(a); hold on; plot(b);
%}

for i=1:size(grt_x,2)
    %--- transition errors
    mt_ekf(i) = norm(ekf_x(1:3,i)-grt_x(1:3,i))./norm(grt_x(1:3,i))*100;
    mt_cepnp(i) = norm(cepnp_x(1:3,i)-grt_x(1:3,i))./norm(grt_x(1:3,i))*100;
    mt_mlpnp(i) = norm(mlpnp_x(1:3,i)-grt_x(1:3,i))./norm(grt_x(1:3,i))*100;
    
    %--- rotatin error
    ekf_dq = quatmultiply(quatinv(grt_x(4:7,i)'), ekf_x(4:7,i)');
    [ekf_dq_yaw, ekf_dq_pitch, ekf_dq_roll] = quat2angle(ekf_dq);
    mr_ekf(i) = norm([ekf_dq_yaw, ekf_dq_pitch, ekf_dq_roll]);
    
    cepnp_dq = quatmultiply(quatinv(grt_x(4:7,i)'), cepnp_x(4:7,i)');
    [cepnp_dq_yaw, cepnp_dq_pitch, cepnp_dq_roll] = quat2angle(cepnp_dq);
    mr_cepnp(i) = norm([cepnp_dq_yaw, cepnp_dq_pitch, cepnp_dq_roll]);

    mlpnp_dq = quatmultiply(quatinv(grt_x(4:7,i)'), mlpnp_x(4:7,i)');
    [mlpnp_dq_yaw, mlpnp_dq_pitch, mlpnp_dq_roll] = quat2angle(mlpnp_dq);
    mr_mlpnp(i) = norm([mlpnp_dq_yaw, mlpnp_dq_pitch, mlpnp_dq_roll]);

end
ekf_err=[mt_ekf; mr_ekf];
cepnp_err=[mt_cepnp; mr_cepnp];
mlpnp_err=[mt_mlpnp; mr_mlpnp];
end
%%
function plotErros(ekf_err, cepnp_err, mlpnp_err)
figure;
subplot(2,1,1);
hold on;
grid on;
plot(ekf_err(1,:), '-r', 'linewidth', 1);
%plot(cepnp_err(1,:), '-b', 'linewidth', 1);
plot(mlpnp_err(1,:), '-k', 'linewidth', 1);
title('Mean Translation Error');
xlabel('n');
ylabel('Translation Error (%)');
legend({'EKFPnP', 'MLPnP'}, 'Location','northwest');
%legend('EKFPnP', 'CEPnP', 'MLPnP');

subplot(2,1,2);
hold on;
grid on;
plot(ekf_err(2,:), '-r', 'linewidth', 1);
%plot(cepnp_err(2,:), '-b', 'linewidth', 1);
plot(mlpnp_err(2,:), '-k', 'linewidth', 1);
title('Mean Rotation Error');
xlabel('n');
ylabel('Rotation Error (degrees)');
legend({'EKFPnP', 'MLPnP'}, 'Location','northwest');
%legend('EKFPnP', 'CEPnP', 'MLPnP');
end