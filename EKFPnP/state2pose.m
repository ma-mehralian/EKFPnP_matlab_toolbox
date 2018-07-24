function P=state2pose(S)
% Conver state vector to Pose
% T = C
% R = R'(q)
%
% INPUT
%   S[7,n]
%
% OUTPUT
%   P[3,4,n]
%
% AUTHOR: ma.mehralian

n=size(S,2);
P=zeros(3,4,n);
for i=1:n
    R=quat2rotm(quatnormalize(S(4:7,i)')); %quat2rotm(S(4:7, i)')';
    t=S(1:3, i);
    P(:,:,i)=[R t];
end
end
