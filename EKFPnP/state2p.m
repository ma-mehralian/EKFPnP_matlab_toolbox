function P=state2p(S)
% Conver state vector to Projection matrix
% P=[R -RC]
%
% INPUT
%   S[7,n]
%
% OUTPUT
%   P[3,4,n]
%
% AUTHOR: ma.mehralian
% Ref: MVG pp156

n=size(S,2);
P=zeros(3,4,n);
%---- make scaler part positive before coversion
si=repmat(sign(S(4,:)), 4, 1);
P(:,1:3,:)=quat2rotm(quatnormalize(S(4:7,:)')); %quat2rotm((si.*S(4:7,:))');
for i=1:n
    P(:,4,i)=-P(:,1:3,i)*S(1:3,i);
end
end
