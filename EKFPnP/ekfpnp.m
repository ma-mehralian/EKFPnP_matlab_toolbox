classdef ekfpnp
    %EKFPnP Extended Kaman Filter for camera pose estimation
    % 
    % .====================== NOTE ======================.
    % |       Some part of this code is based on         |
    % | Javier Civera and J. M. M. Montiel great work on |
    % |        1-Point RANSAC for EKF Filtering          |
    % '=================================================='
    
    properties
        s_k_k
        s_k_km1
        p_k_k
        p_k_km1
        std_a
        std_alpha
        dist
    end
    
    methods
        function obj = ekfpnp(s1, s2, dist, varargin)
            % s1, s2 is the first 2 states of system
            % dist = [f cx, cy]
            % varargin = frames distnce of first 2 states
            
            if ~isempty(varargin)
                ratio = varargin{1};
            else
                ratio = 1;
            end
            % --- initial linear velocity
            v_0 =(s2(1:3)-s1(1:3))/ratio;
            std_v_0 = 0.025;
            
            r=quatmultiply(quatinv(s1(4:7)'), s2(4:7)')';
            theta = 2*acos(r(1));
            if norm(r(2:4)) ~=0
                w_0 = theta*r(2:4)/norm(r(2:4)) /ratio;
            else
                w_0 = [eps eps eps]';
            end
            std_w_0 = 0.025;
            %-----
            %[s1(1:3)+v_0 ; quatmultiply(s1(4:7)',v2q(w_0))'] - s2
            
            % ---- Initialize state vector and covariance
            s_k_k = [s2; v_0; w_0];
            p_k_k = zeros(13,13);
            p_k_k(1:14:7*13)=eps;
            p_k_k(7*13+8:14:10*13)=std_v_0^2;
            p_k_k(10*13+11:14:13*13)=std_w_0^2;
            
            % ---- create EKF object
            obj.s_k_k = s_k_k;
            obj.p_k_k = p_k_k;
            obj.std_a = 0.05;
            obj.std_alpha = 0.05;
            obj.dist = dist;
        end
        %%
        function obj = step_iterated(obj, X, z, varargin)
            if nargin>3
                R = varargin{1};
            else
                R = eye(length(z));
            end
            [ obj.s_k_km1, obj.p_k_km1 ] = ekf_predict( obj.s_k_k, obj.p_k_k, [obj.std_a, obj.std_alpha]);
            [ obj.s_k_k, obj.p_k_k, K ] = ekf_update_iterated( obj.s_k_km1, obj.p_k_km1, z, X, obj.dist, R);
        end
        %%
        function obj = step_simple(obj, X, z, varargin)
            if nargin>3
                R = varargin{1};
            else
                R = eye(length(z));
            end
            [ obj.s_k_km1, obj.p_k_km1 ] = ekf_predict( obj.s_k_k, obj.p_k_k, [obj.std_a, obj.std_alpha]);
            [ obj.s_k_k, obj.p_k_k, K ] = ekf_update_simple_fast( obj.s_k_km1, obj.p_k_km1, z, X, obj.dist, R);
        end
        %%
        function obj = step_simple_old_slow(obj, X, z, varargin)
            if nargin>3
                R = varargin{1};
            else
                R = eye(length(z));
            end
            [ obj.s_k_km1, obj.p_k_km1 ] = ekf_predict( obj.s_k_k, obj.p_k_k, [obj.std_a, obj.std_alpha]);
            [ obj.s_k_k, obj.p_k_k, K ] = ekf_update_simple_old_slow( obj.s_k_km1, obj.p_k_km1, z, X, obj.dist, R);
        end
        %%
        function [R, t] = getProj(obj)
            P=state2p(obj.s_k_k);
            R=P(:,1:3);
            t=P(:,4);
        end
        %%
        function [R, t] = getPose(obj)
            R=quat2rotm(obj.s_k_k(4:7)')';
            t=obj.s_k_k(1:3);
        end
    end
end

%% ##################### PREDICTION ######################
%%
function [s_k_km1, P_k_km1]=ekf_predict(s_km1_km1, P_km1_km1, u_k)
s_k_km1 = f_s(s_km1_km1, u_k);
% figure; structure_plot(cat(3, state2p(s_k_km1(1:7)),state2p(s_km1_km1(1:7))), [1,1,1]');
F = dfs_by_ds(s_k_km1);
G = dfs_by_du(s_k_km1);

C_a = u_k(1)^2;
C_alpha = u_k(2)^2;
Q=diag([C_a C_a C_a C_alpha C_alpha C_alpha]);
P_k_km1 = F*P_km1_km1*F'+G*Q*G';
end
%%
function s_k_km1=f_s(s_km1_km1, u_k)
dt=1;
C = s_km1_km1(1:3,1);
q = s_km1_km1(4:7,1);
v = s_km1_km1(8:10,1);
V = 0;%randn*u_k(1)*dt;
w = s_km1_km1(11:13,1);
W = 0;%randn*u_k(2)*dt;

s_k_km1=[C+(v+V)*dt;
    quatmultiply(q',v2q((w+W)*dt))';
    (v+V);
    (w+W)];
end
%%
function q=v2q(v)
%converts rotation vector to quaternion.
nrm = norm(v);
if (nrm <eps)
    q=[1 0 0 0];
else
    v_n=reshape(v,[1,3])/nrm;
    theta=nrm;
    q=[cos(theta/2)  sin(theta/2)*v_n];
end
end
%%
function G=dfs_by_du(s_k_km1)
dt=1;
C =s_k_km1(1:3);
q =s_k_km1(4:7);
v =s_k_km1(8:10);
w =s_k_km1(11:13);

G=zeros(13,6);
G(8:10,1:3)=eye(3);
G(11:13,4:6)=eye(3);
G(1:3,1:3)=eye(3)*dt;
G(4:7,4:6)=dq3_by_dq1(q)*dqomegadt_by_domega(w,dt);
end
%%
function F=dfs_by_ds(s_k_km1)
dt=1;
C =s_k_km1(1:3);
q =s_k_km1(4:7);
v =s_k_km1(8:10);
w =s_k_km1(11:13);

F=eye(13);
F(1:3,8:10)= eye(3)*dt;

qwt=v2q(w*dt);
F(4:7,4:7) = dq3_by_dq2(qwt);
F(4:7,11:13) = dq3_by_dq1(q)*dqomegadt_by_domega(w,dt);
end
%%
function dq3_by_dq2RES=dq3_by_dq2(q1_in)
q1.R=q1_in(1);
q1.X=q1_in(2);
q1.Y=q1_in(3);
q1.Z=q1_in(4);

dq3_by_dq2RES = ...
   [q1.R, -q1.X, -q1.Y, -q1.Z,
    q1.X,  q1.R,  q1.Z, -q1.Y,
    q1.Y, -q1.Z,  q1.R,  q1.X,
    q1.Z,  q1.Y, -q1.X,  q1.R];
end
%%
function dq3_by_dq1RES=dq3_by_dq1(q2_in)
q2.R=q2_in(1);
q2.X=q2_in(2);
q2.Y=q2_in(3);
q2.Z=q2_in(4);

dq3_by_dq1RES=[q2.R, -q2.X, -q2.Y, -q2.Z,
    q2.X,  q2.R, -q2.Z,  q2.Y,
    q2.Y,  q2.Z,  q2.R, -q2.X,
    q2.Z, -q2.Y,  q2.X,  q2.R];
end
%%
function dqomegadt_by_domegaRES=dqomegadt_by_domega(omega, dt)

%// Modulus
omegamod = norm(omega);

%// Use generic ancillary functions to calculate components of Jacobian
dqomegadt_by_domegaRES(1, 1) = dq0_by_domegaA(omega(1), omegamod, dt);
dqomegadt_by_domegaRES(1, 2) = dq0_by_domegaA(omega(2), omegamod, dt);
dqomegadt_by_domegaRES(1, 3) = dq0_by_domegaA(omega(3), omegamod, dt);
dqomegadt_by_domegaRES(2, 1) = dqA_by_domegaA(omega(1), omegamod, dt);
dqomegadt_by_domegaRES(2, 2) = dqA_by_domegaB(omega(1), omega(2), omegamod, dt);
dqomegadt_by_domegaRES(2, 3) = dqA_by_domegaB(omega(1), omega(3), omegamod, dt);
dqomegadt_by_domegaRES(3, 1) = dqA_by_domegaB(omega(2), omega(1), omegamod, dt);
dqomegadt_by_domegaRES(3, 2) = dqA_by_domegaA(omega(2), omegamod, dt);
dqomegadt_by_domegaRES(3, 3) = dqA_by_domegaB(omega(2), omega(3), omegamod, dt);
dqomegadt_by_domegaRES(4, 1) = dqA_by_domegaB(omega(3), omega(1), omegamod, dt);
dqomegadt_by_domegaRES(4, 2) = dqA_by_domegaB(omega(3), omega(2), omegamod, dt);
dqomegadt_by_domegaRES(4, 3) = dqA_by_domegaA(omega(3), omegamod, dt);

end


% // Ancillary functions: calculate parts of Jacobian dq_by_domega
% // which are repeatable due to symmetry.
% // Here omegaA is one of omegax, omegay, omegaz
% // omegaB, omegaC are the other two
% // And similarly with qA, qB, qC

function dq0_by_domegaARES=dq0_by_domegaA(omegaA, omega, dt)
dq0_by_domegaARES=(-dt / 2.0) * (omegaA / omega) * sin(omega * dt / 2.0);
end

function dqA_by_domegaARES=dqA_by_domegaA(omegaA, omega, dt)
dqA_by_domegaARES=(dt / 2.0) * omegaA * omegaA / (omega * omega) ...
    * cos(omega * dt / 2.0) ...
    + (1.0 / omega) * (1.0 - omegaA * omegaA / (omega * omega))...
    * sin(omega * dt / 2.0);
end

function dqA_by_domegaBRES=dqA_by_domegaB(omegaA, omegaB, omega, dt)
dqA_by_domegaBRES=(omegaA * omegaB / (omega * omega)) * ...
    ( (dt / 2.0) * cos(omega * dt / 2.0) ...
    - (1.0 / omega) * sin(omega * dt / 2.0) );
end
%% ##################### CORRECTION ######################
%%
function [ s_k_k, p_k_k, K ] = ekf_update_iterated( s_km1_k, p_km1_k, z, X, dist, R)
s_k_k=s_km1_k;
for i=1:2
    [h, X_C]=f_h(s_k_k, X, dist);
    H = dfh_by_ds(s_k_k, X, X_C, dist);
    if i==1, h1=h; end
    
    % --- filter gain
    S = H*p_km1_k*H' + R;
    K = p_km1_k*H'/S;
    
    % --- updated state
    s_k_k = s_km1_k + K*( z - h -H*(s_km1_k-s_k_k) );
end
% --- updated covariance
p_k_k = (eye(13) - K*H)*p_km1_k;

% --- normalize the quaternion
Jnorm = normJac( s_k_k( 4:7 ) );
size_p_k_k = size(p_k_k,1);
p_k_k = [   p_k_k(1:3,1:3)              p_k_k(1:3,4:7)*Jnorm'               p_k_k(1:3,8:size_p_k_k);
    Jnorm*p_k_k(4:7,1:3)        Jnorm*p_k_k(4:7,4:7)*Jnorm'         Jnorm*p_k_k(4:7,8:size_p_k_k);
    p_k_k(8:size_p_k_k,1:3)     p_k_k(8:size_p_k_k,4:7)*Jnorm'      p_k_k(8:size_p_k_k,8:size_p_k_k)];
s_k_k( 4:7 ) = s_k_k( 4:7 ) / norm( s_k_k( 4:7 ) );


end
%%
function [ s_k_k, p_k_k, K ] = ekf_update_simple_fast( s_km1_k, p_km1_k, z, X, dist, R)
[h, X_C]=f_h(s_km1_k, X, dist);
H = dfh_by_ds(s_km1_k, X, X_C, dist);

% --- filter gain
inv_R_D = diag(R).^-1;
%B = (p_km1_k * H')* diag(inv_R_D);
B = bsxfun(@times,(p_km1_k * H'), inv_R_D');
C = B*H;
K = B - C *inv(eye(13) + C )*B;

% --- updated state and covariance
s_k_k = s_km1_k + K*( z - h );
p_k_k = (eye(13) - K*H)*p_km1_k;

% --- normalize the quaternion
Jn = normJac( s_k_k( 4:7 ) );
p_k_k = [...
    p_k_k(1:3,1:3)      zeros(3,4)       p_k_k(1:3,8:10)         zeros(3);
    zeros(4,3)   Jn*p_k_k(4:7,4:7)*Jn'     zeros(4,3)       Jn*p_k_k(4:7,11:13);
    p_k_k(8:10,1:3)           zeros(3,4)        p_k_k(8:10,8:10)        zeros(3);
    zeros(3)      p_k_k(11:13,4:7)*Jn'           zeros(3)         p_k_k(11:13,11:13)];
s_k_k( 4:7 ) = s_k_k( 4:7 ) / norm( s_k_k( 4:7 ) );
end
%%
function [ s_k_k, p_k_k, K ] = ekf_update_simple_old_slow( s_km1_k, p_km1_k, z, X, dist, R)
[h, X_C]=f_h(s_km1_k, X, dist);
H = dfh_by_ds(s_km1_k, X, X_C, dist);

% --- filter gain
S = H*p_km1_k*H' + R;
K = p_km1_k*H'/S;

% plot(z(1:2:end), z(2:2:end), '.b'); hold on; plot(h(1:2:end), h(2:2:end), '.r');
% --- updated state and covariance
s_k_k = s_km1_k + K*( z - h );
p_k_k = (eye(13) - K*H)*p_km1_k;

% --- normalize the quaternion
Jnorm = normJac( s_k_k( 4:7 ) );
size_p_k_k = size(p_k_k,1);
p_k_k = [   p_k_k(1:3,1:3)              p_k_k(1:3,4:7)*Jnorm'               p_k_k(1:3,8:size_p_k_k);
    Jnorm*p_k_k(4:7,1:3)        Jnorm*p_k_k(4:7,4:7)*Jnorm'         Jnorm*p_k_k(4:7,8:size_p_k_k);
    p_k_k(8:size_p_k_k,1:3)     p_k_k(8:size_p_k_k,4:7)*Jnorm'      p_k_k(8:size_p_k_k,8:size_p_k_k)];
s_k_k( 4:7 ) = s_k_k( 4:7 ) / norm( s_k_k( 4:7 ) );


end
%%
function [h, X_C]=f_h(s, X, dist)
[x, X_C]=qproj(s, dist, X);
h=reshape(x, [], 1);
end
%%
function H=dfh_by_ds(s, X, X_C, dist)
% Projection function: h= h_CI(h_WC) : world -> camera -> image
% dh/ds = d(h_CI)/d(H_WC) * d(H_WC)/ds
% ---
C=s(1:3);
q=s(4:7);
qr=q(1);
qv=q(2:4);
n=size(X_C,2);
H = zeros(2*n, 13);
f=dist(1);

% --- d(h_CI)/d(H_WC)
dhci_by_dhwc = zeros(2*n,3);
dhci_by_dhwc(1:2:end,:) =f*[1./X_C(3,:)' zeros(n,1) -(X_C(1,:)./(X_C(3,:).^2))'];
dhci_by_dhwc(2:2:end,:) =f*[zeros(n,1) 1./X_C(3,:)' -(X_C(2,:)./(X_C(3,:).^2))'];

% --- d(H_WC)/ds
dhwc_by_dc = -( (qr^2-norm(qv)^2)*eye(3,3) + 2*qv*qv' + 2*qr*skew(qv) );
% dhwc_by_dc = -( quat2rotm(q) );

% ---
v=X-repmat(C, 1, size(X,2));
qv_rep=repmat(qv,[1,size(X,2)]);
dot_v_qv=sum(qv_rep.*v);
crs_v_qv=qr*v+cross(qv_rep, v, 1);

for i=1:n
    vi= v(:,i);
    %{
    dwc_by_dq = 2*[q0*vi + cross(qv,vi), qv'*vi*eye(3,3) + qv*vi' - vi*qv' - q0*skew(vi)];
    %}
    tmp=qv*vi';
    q0_pv=qr*vi;
    dwc_by_dq = 2*[crs_v_qv(:,i), (dot_v_qv(i)*eye(3,3)+ tmp -tmp' - [0 -q0_pv(3) q0_pv(2); q0_pv(3) 0 -q0_pv(1); -q0_pv(2) q0_pv(1) 0]) ];
    %---
    
    dhh_by_dx_m=[dhwc_by_dc dwc_by_dq zeros(3,6)];
    ix=2*(i-1)+1;
    dhci_by_dhwc_m = dhci_by_dhwc(ix:ix+1,:);
    H(ix:ix+1, :) = dhci_by_dhwc_m*dhh_by_dx_m;
end
end
%%
function M = skew(t)
M = [0 -t(3) t(2)
    t(3) 0 -t(1)
    -t(2) t(1) 0];
end
%%
function J=normJac(q)
r=q(1);
x=q(2);
y=q(3);
z=q(4);

J=(r*r+x*x+y*y+z*z)^(-3/2)*...
    [x*x+y*y+z*z         -r*x         -r*y         -r*z;
    -x*r  r*r+y*y+z*z         -x*y         -x*z;
    -y*r         -y*x  r*r+x*x+z*z         -y*z;
    -z*r         -z*x         -z*y  r*r+x*x+y*y];
end