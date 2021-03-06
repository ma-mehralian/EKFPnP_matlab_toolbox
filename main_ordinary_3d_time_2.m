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

clear; clc;
IniToolbox;

disp(['START: ', datestr(datetime('now'))]);
warning('off');

% experimental parameters
nl= 1:10;
nlsamples = 0.1*ones(1,10); %percentatge of samples for each sigma
npts= [10,100:100:2000];
num = 50;
cam_num = 30;

% compared methods
A= zeros(size(npts));
B= zeros(num,1);

name = {'slow-EKFPnP', 'EKFPnP', 'MLPnP','EPnP+GN', 'RPnP', 'ASPnP', 'OPnP', 'DLS', 'CEPPnP', 'SRPnP'};
f = {[], [], @MLPNP_with_COV, @EPnP_GN,  @RPnP, @ASPnP, @OPnP, @robust_dls_pnp, @CEPPnP, @SRPnP};
marker = { 'o', 'x', 'd', '*', 's', 'v', 'o', '+', '>', '^'};
color = {[.5,0.3,0], [.8,0,0], [0,0,.8], [.6,.4,.2], [.8,0,.8], [0,0,0], [.2,.5,.8], [1,.5,0], [0.8,.5,0.1], [0.1,.5,0.1]};
markerfacecolor = color;


method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor);



meanR = zeros(length(npts),length(method_list)+1);
medianR = zeros(length(npts),length(method_list)+1);
meanT = zeros(length(npts),length(method_list)+1);
medianT = zeros(length(npts),length(method_list)+1);
medianTime = zeros(length(npts),length(method_list)+1);
% experiments
for i= 1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d (num sg = %d ): ', npt, length(nl));
    
    for k= 1:length(method_list)
        method_list(k).c = zeros(1,num);
        method_list(k).e = zeros(1,num);
        method_list(k).r = zeros(1,num);
        method_list(k).t = zeros(1,num);
    end
    index_fail = cell(1,length(name));
    
    for j= 1:num
        % camera's parameters
        width= 640;
        height= 480;
        f= 800;
        K = [f 0 0
            0 f 0
            0 0 1];
        
        % generate camera poses in world coordinate
        Ps=sample_seq('cam_num', 200, 'traj_r', 10);
        Ps=Ps(:,:,1:cam_num);
        % generate 3d points in world coordinate
        XXw= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
        
        nnl = round(npt * nlsamples);
        nls = zeros(1,npt);
        id  = 1;
        for idnl = 1:length(nl)
            sigma = nl(idnl);
            nls(id:id+nnl(idnl)-1) = sigma .* ones(1,nnl(idnl));
            id = id+nnl(idnl);
        end
        
        %--- for each camera in sequense Ps
        for h=1:cam_num
            inx = (j-1)*cam_num+h;
            P = pInv(Ps(:,:,h));
            t = P(:,4);
            R = P(:,1:3);
            Xc= P*[XXw;ones(1,npt)];
            
            % projection
            xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
            
            
            randomvals = randn(2,npt);
            nl_e = log10(1+100*(h-1)/cam_num)/2;
            xxn= xx+randomvals.*(nl_e*[nls;nls]);
            
            homx = [xxn/f; ones(1,size(xxn,2))];
            v = normc(homx);
            
            Cu = zeros(2,2,npt);
            Evv = zeros(3,3,npt);
            cov = zeros(9,size(Cu,3));
            for id = 1:npt
                Cu(:,:,id) = diag([mean(nl_e*nl)^2 mean(nl_e*nl)^2]);
                cov_proj = K\diag([mean(nl_e*nl)^2 mean(nl_e*nl)^2 0])/K';
                J = (eye(3)-(v(:,id)*v(:,id)')/(v(:,id)'*v(:,id)))/norm(homx(:,id));
                Evv(:,:,id) = J*cov_proj*J';
                cov(:,id) = reshape(Evv(:,:,id),9,1);
            end
            RC = mean(nl_e*nl)*eye(2*npt);
            
            
            % pose estimation
            for k= 1:length(method_list)
                try
                    if strcmp(method_list(k).name, 'CEPPnP')
                        tic;
                        mXXw = XXw - repmat(mean(XXw,2),1,size(XXw,2));
                        [R1,t1]= method_list(k).f(mXXw,xxn/f,Cu);
                        t1 = t1 - R1 * mean(XXw,2);
                        tcost = toc;
                    elseif strcmp(method_list(k).name, 'MLPnP') || strcmp(method_list(k).name, 'MLPnPWithCov') 
                        tic;
                        [R1,t1]= method_list(k).f(XXw,v,cov);
                        tcost = toc;
                    elseif strcmp(method_list(k).name, 'EKFPnP')
                        tic;
                        %--- initialize the first two view with EPnP_GN
                        if h<=2
                            [R1,t1] = EPnP_GN(XXw, xxn/f);
                            INT_P(:,:,h)= pInv([R1,t1]);
                            if h==2 , ekf_obj = ekfpnp(pose2state(INT_P(:,:,1)), pose2state(INT_P(:,:,2)), [f 0 0]); end
                        elseif h>2
                            z= reshape(xxn,[],1);
                            ekf_obj = ekf_obj.step_simple(XXw, z, RC);
                            [R1, t1] = ekf_obj.getProj();
                            
                        end
                        tcost = toc;
                    elseif strcmp(method_list(k).name, 'slow-EKFPnP')
                        tic;
                        %--- initialize the first two view with EPnP_GN
                        if h<=2
                            [R1,t1] = EPnP_GN(XXw, xxn/f);
                            slow_INT_P(:,:,h)= pInv([R1,t1]);
                            if h==2 , slow_ekf_obj = ekfpnp(pose2state(slow_INT_P(:,:,1)), pose2state(slow_INT_P(:,:,2)), [f 0 0]); end
                        elseif h>2
                            z= reshape(xxn,[],1);
                            slow_ekf_obj = slow_ekf_obj.step_simple_old_slow(XXw, z, RC);
                            [R1, t1] = slow_ekf_obj.getProj();
                            
                        end
                        tcost = toc;
                    else
                        tic;
                        [R1,t1]= method_list(k).f(XXw,xxn/f);
                        tcost = toc;
                    end
                catch
                    %disp(['The solver - ',method_list(k).name,' - encounters internal errors!!!']);
                    %index_fail = [index_fail, j];
                    index_fail{k} = [index_fail{k}, inx];
                    break;
                end
                
                %no solution
                if size(t1,2) < 1
                    disp(['The solver - ',method_list(k).name,' - returns no solution!!!']);
                    %index_fail = [index_fail, j];
                    index_fail{k} = [index_fail{k}, inx];
                    break;
                elseif (sum(sum(sum(imag(R1).^2))>0) == size(R1,3) || sum(sum(imag(t1(:,:,1)).^2)>0) == size(t1,2))
                    index_fail{k} = [index_fail{k}, inx];
                    continue;
                end
                %choose the solution with smallest error
                error = inf;
                for jjj = 1:size(R1,3)
                    tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                    if sum(tempy) < error
                        cost  = tcost;
                        ercorr= mean(sqrt(sum((R1(:,:,jjj) * XXw +  t1(:,jjj) * ones(1,npt) - Xc).^2)));
                        y     = tempy;
                        error = sum(tempy);
                    end
                end
                
                method_list(k).c(inx)= cost * 1000;
            end
            showpercent(inx,num*cam_num);
        end
    end
    fprintf('\n');
    
    % save result
    for k= 1:length(method_list)
        %results without deleting solutions
        tmethod_list = method_list(k);
        method_list(k).c(index_fail{k}) = [];
        
        % computational cost should be computed in a separated procedure as
        % in main_time.m
        
        method_list(k).pfail(i) = 100 * numel(index_fail{k})/(num*cam_num);
        
        method_list(k).mean_c(i)= mean(method_list(k).c);
        method_list(k).med_c(i)= median(method_list(k).c);
        method_list(k).std_c(i)= std(method_list(k).c);


        medianTime(i,1) = npts(i);      
        medianTime(i,k+1) = method_list(k).med_c(i);
        
        %results deleting solutions where not all the methods finds one
        tmethod_list.c(unique([index_fail{:}])) = [];
        
        method_list(k).deleted_mean_c(i)= mean(tmethod_list.c);
        method_list(k).deleted_med_c(i)= median(tmethod_list.c);
        method_list(k).deleted_std_c(i)= std(tmethod_list.c);       
    end
end

save ordinary3DresultsTime method_list npts;
disp(['END: ', datestr(datetime('now'))]);
