function P_inv=pInv(P)
% Inverse Pose/Projection Matrix (projection2pose and pose2projection)
% INPUT:
%   P[3x4xn]
% OUTPUT
%   P_inv[3x4xn]

n=size(P,3);
P_inv=zeros(3,4,n);
for i=1:n
    R=P(:,1:3,i);
    t=P(:,4,i);
    P_inv(:,:,i)=[R' -R'*t];
end
end