function [Ps, Ps_q]=sample_seq(varargin)
% Genarate circular camera trajectory
%
% OUTPUT:
%   Ps[3,4,n]: n camera pose ([R,t])
%   Ps_q[7,n]: n camera pose ([t;q(R)])
%
% OPTIONS:
%   'cam_num', The number of cameras. Default= 150
%   'traj_r', camera trajectory radius. Default= 10
%
% AUTHOR: ma.mehralian

opts=pars_options(varargin{:});

r=opts.traj_r;
n=opts.cam_num;
ang=linspace(0, 6*pi, n)';
Ps=zeros(3,4,n);
Ps_q=zeros(7,n);
for i=1:n
    %--- g=[-0.5,0.5]
    g= i/n-0.5;
    r_ = r*(1-abs(1.3*g));
    center=[r_*sin(ang(i)) r_*cos(ang(i)) 15*g+6]';
    a = -ang(i);%0;
    b = pi/2 + g*pi/2;%pi/2;
    c = -ang(i);%0;
    
    q1=[cos(a/2) 0 0 sin(a/2)];
    q2=[cos(b/2) sin(b/2) 0 0];
    q3=[cos(c/2) 0 0 sin(c/2)];
    q=quatmultiply(quatmultiply(q1,q2) ,q3);
    
    Ps(:,:,i)=[quat2rotm(q) center];
    Ps_q(:,i)=[center; q'];
end
% npt = 300;
% plot_cam_pose(Ps, [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])], 'cam_scale', 0.3);
end
%%
function opts=pars_options(varargin)
%--- defaults
opts.cam_num    = 150;
opts.traj_r     = 10;
%--- start parse
if nargin > 0
    i = 1;
    while i <= nargin
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'cam_num'
                    opts.cam_num = varargin{i+1};
                    i=i+1;
                case 'traj_r'
                    opts.traj_r = varargin{i+1};
                    i=i+1;
            end
        end
        i = i + 1;
    end
end
end