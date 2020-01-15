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

function main_box
clear; clc;
IniToolbox;

id = 'b3';
video = VideoReader(['data/', id ,'.mp4']);
% vo = VideoWriter('box_results.avi');
% open(vo);

load(['data/',id,'.mat']);
pts_2d = pts_2d(:,1:20:end);
pts_3d = pts_3d(:,1:20:end);

invK= inv(K);
f = K(1,1);
c = K(1:2,3);
[R1, t1] = efficient_pnp_gauss(pts_3d', pts_2d', K);
edges=[vs(:,1) vs(:,2) vs(:,6) vs(:,5) vs(:,1) ...
    vs(:,3) vs(:,7) vs(:,5) vs(:,6) vs(:,8) vs(:,4) ...
    vs(:,3) vs(:,7) vs(:,8) vs(:,6) vs(:,2) vs(:,4)];

img1 = im2single(rgb2gray(readFrame(video)));
tracker = vision.PointTracker();
initialize(tracker, pts_2d', img1);
ekf_initialized = false;
i=0;
figure('outerposition',[0  100  1920  400])
while hasFrame(video)
    i=i+1;
    disp(sprintf('frame %d', i));
    frame = readFrame(video);
    imgi= im2single(rgb2gray(frame));
    [xi, validity] = step(tracker, imgi);
    
    X=pts_3d(:,validity);
    x=xi(validity,:)';
    npt = size(x,2);
    x_i = (x-repmat(c,1,npt))/f;
    homx = [x_i; ones(1,npt)];
    v = normc(homx);
    
    %{
    %--- EPnP+GN/RPnP/OPnP
    [R_epnp,t_epnp] = EPnP_GN(X, x_i);
    [R_opnp,t_opnp] = OPnP(X, x_i);
    [R_rpnp,t_rpnp] = robust_dls_pnp(X, x_i);
    %}
    
    Cu = zeros(2,2,npt);
    Evv = zeros(3,3,npt);
    cov = zeros(9,size(Cu,3));
    nl = 3;
    for id = 1:npt
        Cu(:,:,id) = diag([nl^2 nl^2]);
        cov_proj = K\diag([nl^2 nl^2 0])/K';
        J = (eye(3)-(v(:,id)*v(:,id)')/(v(:,id)'*v(:,id)))/norm(homx(:,id));
        Evv(:,:,id) = J*cov_proj*J';
        cov(:,id) = reshape(Evv(:,:,id),9,1);
    end
    %--- CEPNP
    %Cu = repmat(3^2 *eye(2),[1,1,size(X,2)]);
    mX = X - repmat(mean(X,2),1,size(X,2));
    [R_cepnp,t_cepnp]=CEPPnP(mX, x_i, Cu);
    t_cepnp = t_cepnp - R_cepnp * mean(X,2);
    
    %--- MLPnP
    v = normc([x_i; ones(1,size(x_i,2))]);
    [R_mlpnp,t_mlpnp]=MLPNP_without_COV(X, v, []);

    
    %--- EKFPnP
    if ~ekf_initialized
        [R_epnp,t_epnp] = EPnP_GN(X, x_i);
        P1 = [R1, t1];
        P2 = [R_epnp,t_epnp];
        ekf_obj = ekfp(p2state(P1), p2state(P2), [K(1,1) K(1:2,3)']);
        ekf_initialized = true;
        R_ekf=R_epnp; t_ekf=t_epnp;
    else
        z=reshape(x,[],1);
        R = 3^2 *eye(length(z));
        ekf_obj = ekf_obj.step_simple(X, z);
        [R_ekf, t_ekf] = ekf_obj.getProj();
    end

    %------- plot results
    %{
    subplot(2,3,1)
    plot_box(K, R_epnp, t_epnp, x, edges, imgi);
    title('EPnP+GN');

    
    subplot(2,3,2)
    plot_box(K, R_rpnp, t_rpnp, x, edges, imgi);
    title('DLS');
    
    subplot(2,3,3)
    plot_box(K, R_opnp, t_opnp, x, edges, imgi);
    title('OPnP');
    %}
    
    subplot(1,9,[1:3])
    plot_box(K, R_cepnp, t_cepnp, x, edges, imgi, i);
    title('CEPnP');

    subplot(1,9,[4:6])
    plot_box(K, R_mlpnp, t_mlpnp, x, edges, imgi, i);
    title('MLPnP');
    
    subplot(1,9,[7:9])
    plot_box(K, R_ekf, t_ekf, x, edges, imgi, i);
    title('EKFPnP');
    
    drawnow;
    %F = getframe(gcf);
    %writeVideo(vo,F.cdata);
    
%     plot_box(K, R_ekf, t_ekf, x, edges, imgi);
%     set(gcf, 'outerposition',[400  200  500  360])
%     drawnow;
    if ismember(i,[1, 400, 450, 500, 550])
        saveas(gcf, sprintf('real_experiment_%d', i) ,'epsc')
    end
end
% close(vo);
end
%%
function plot_box(K, R, t, x, e, frame, frame_index)
% K:calibration matrix
% R,t : camera pose
% x: 2D Features
% e: box edges
    e_new = K*[R t]* [e; ones(1,length(e))];
    e_new = [e_new(1,:)./e_new(3,:); e_new(2,:)./e_new(3,:)];
    hold off; imshow(frame); hold on; 
    plot(x(1,:), x(2,:),'.r');
    %plot(e_new(1,:), e_new(2,:),'.r');
    line(e_new(1,:), e_new(2,:), 'LineWidth', 2, 'Color', [0 1 0]);
    text(20,40,sprintf('Frame %d', frame_index), 'FontSize',12)
end