function S=p2state(P)
% Conver Projection matrix to state vector
% C= -R'T
% R= R
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
    Rq=rotm2quat(P(:, 1:3, i));
    %---- make normal and scaler part positive
    Rq=Rq/norm(Rq);% * sign(Rq(1));
    
    C=-P(:, 1:3, i)'*P(:,4,i);
    S(:,i)=[C; Rq'];
end
end
