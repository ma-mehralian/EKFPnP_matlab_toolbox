function S=pose2state(P)
% Conver Pose to state vector
% C = t
% q = q(R')
%
% INPUT
%   P[3,4,n]
%
% OUTPUT
%   S[7,n]
%
% AUTHOR: ma.mehralian

n=size(P,3);
S=zeros(7,n);
for i=1:n
    q=quatnormalize(rotm2quat(P(:, 1:3, i)')); %rotm2quat(P(:, 1:3, i)');   
    C=P(:,4,i);
    S(:,i)=[C; q'];
end
end
