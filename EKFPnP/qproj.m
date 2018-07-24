function [x, X_C] = qproj(s, k, X_W)
% world to image projection (3D to 2D) using camera state in quaternion representation 
%
% INPUT:
%   K[3]  : camera parameters [f, c_i, c_j]
%   s[7]  : camera state
%   X_W[3,n]: n 3D points in world coordinate
%
% OUTPUT:
%   x[2,n,m] : n 2D points for m camera in image coordinate
%   X_C[3,n,m]: n 3D points for m camera in camera coordinate
%              NOTE: this output used in ekf_update

C=s(1:3);
q=s(4:7)';
n=size(X_W,2);
f=k(1);
c=k(2:3)';

v = zeros(4,n);
v(2:end,:) = X_W-repmat(C,1,n);

% X_C_1=quatrotate(quatconj(q), v(2:4,:)')'; %?????
% X_C_2=quat2rotm(q)*v(2:4,:); %?????
X_C=quatmultiply(quatmultiply(q, v'), quatinv(q))';
X_C=X_C(2:4,:);
 
x=f*X_C(1:2,:) ./ repmat(X_C(3,:),[2,1]);
x=repmat(c, [1,n])+x;
end