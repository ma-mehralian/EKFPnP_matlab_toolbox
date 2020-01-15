function [F, inliers, e]= fundamental_projective(x1, x2, varargin)
% FUNDAMENTAL: calculate fundamental matrix
%
% INPUTS:
%   x1[3,m]: homogenous points in first image
%   x2[3,m]: correspondence homogenous points in second image
%
% OPTIONS:
%   'is_calibrated', indicates the camera is calibrated.
%   'use_ransac', whether to use RANSAC for outlier rejection or not.
%   'max_iter', maximum iteration for RANSAC [1000]
%   'dis_thresh', Distance threshold [0.0005]
%   'stop_per', minimum percent of inliers to stop RANSAC [0.95]
%
% OUTPUT:
%   F[3,3]: Fundamental matrix
%   inliers: RANSAC inliers
%   3[3,1]: Right null-vector of F
%
% REFERENCE:
%   HZ2, Algorithm 11.1. The normalized 8-point algorithm for F, p 282
%   HZ2, Algorithm 4.2. The normalized DLT, Page 109
%
% AUTHOR: ma.mehralian

opts=pars_options(varargin{:});

%--- normalize points
[nx1, T1] = normalise2dpts(x1);
[nx2, T2] = normalise2dpts(x2);

if opts.use_ransac
    [F, inliers]=CalFundamentalRANSAC(nx1, nx2, opts);
else
    inliers = 1:length(nx1);
    F = CalFundamental(nx1, nx2, opts);
end

%--- Denormalized fundamental matrix
F=T2'*F*T1;

e= null(F');
%F=F./F(3,3);
end
%%
function F=CalFundamental(nx1, nx2, opts)
%--- Build the constraint matrix (Kronecker product)
A = [nx2(1,:)'.*nx1(1,:)'	nx2(1,:)'.*nx1(2,:)'	nx2(1,:)'.*nx1(3,:)' ...
    nx2(2,:)'.*nx1(1,:)'	nx2(2,:)'.*nx1(2,:)'	nx2(2,:)'.*nx1(3,:)' ...
    nx2(3,:)'.*nx1(1,:)'	nx2(3,:)'.*nx1(2,:)'	nx2(3,:)'.*nx1(3,:)'];
[~,~,V] = svd(A,0);
F_ = reshape(V(:,end),3,3)';

[U,D,V] = svd(F_,0);
%--- Enforce rank 2 constraint and s1 = s2 constraint (for calibrated)
if opts.is_calibrated
    s = (D(1,1)+D(2,2))/2;
    D = diag([s s 0]);
else
    D = diag([D(1,1), D(2,2) 0]);
end
F = U*D*V';
end
%%
function [F, inliers]=CalFundamentalRANSAC(nx1, nx2, opts)
inliers=[];
for i=1:opts.max_iter
    % ---- STEP1: Select at random a set of 8 matches.
    ixs= randsample(length(nx1), 8);
    p1= nx1(:,ixs);
    p2= nx2(:,ixs);
    
    % ---- STEP2: Compute the fundamental matrix F
    c_F = CalFundamental(p1, p2, opts);
    
    % ---- STEP3: Determine the subset of inliers using Sampson distance
    dis = Sampson(nx1, nx2, c_F);
    c_inliers = dis < opts.dis_thresh;
    
    % ---- STEP4: Count the number of points in the consensus set.
    if (sum(c_inliers)>sum(inliers))
        inliers= c_inliers;
        F= c_F;
        if (opts.stop_per <= sum(inliers)/length(inliers) )
            break;
        end
    end
end
F = CalFundamental(nx1(:,inliers), nx2(:,inliers), opts);
end
%%
function d = Sampson(x1, x2, F)
% SAMPSON calculate sampson distance
%   x1: 3xm homogenous points in first image
%   x2: 3xm homogenous points in first image
%   F : 3x3 fundamental matrix calculated from points x1,x2

x2tFx1 = zeros(1,length(x1));
for n = 1:length(x1)
    x2tFx1(n) = x2(:,n)'*F*x1(:,n);
end

Fx1 = F*x1;
Ftx2 = F'*x2;

% Evaluate distances
d =  x2tFx1.^2 ./ ...
    (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);
end
%%
function opts=pars_options(varargin)
%--- defaults
opts.is_calibrated = 0;
opts.use_ransac = 0;
opts.max_iter = 1000;
opts.dis_thresh = 5e-4;
opts.stop_per = 0.95;
%--- start parse
if nargin > 0
    i = 1;
    while i <= nargin
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'is_calibrated'
                    opts.is_calibrated = 1;
                case 'use_ransac'
                    opts.use_ransac = 1;
                case 'max_iter'
                    opts.max_iter = varargin{i+1};
                    i=i+1;
                case 'dis_thresh'
                    opts.dis_thresh = varargin{i+1};
                    i=i+1;
                case 'stop_per'
                    opts.stop_per = varargin{i+1};
                    i=i+1;
            end
        end
        i = i + 1;
    end
end
end