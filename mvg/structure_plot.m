function structure_plot(P, X, varargin)
% Plot reconstructed structure
%
% INPUT:
%   P[3,4,m]: m camera matrices
%   X[3,n]  : n 3D points
%   varargin
%       - 'k'[3,3]  : calibration matrix
%       - 'imgs'[h,w,n]: m camera images
%       - 'cam_color': camera color
%       - 'cam_fill': camera fill color
%       - 'ply_file': ply export of structure plot
%
% AUTHOR: ma.mehralian

opts=pars_options(varargin{:});
n=size(X,2);
m=size(P,3);
sz=1;%sign(mean(X(3,:)));
if size(X,1)==4
    X=[X(1,:)./X(4,:); X(2,:)./X(4,:); sz*X(3,:)./X(4,:)];
end
%scatter3(X(1,:), X(2,:), X(3,:),'filled');
plot3(X(1,:), X(2,:), X(3,:),'.');
grid on;

faces=[];
edges=[];
vertices=[X(1:3,:)' zeros(n,3)];
%return
hold all;
pts = [0 0 0; 1 1 sz; 1 -1 sz; -1 -1 sz; -1 1 sz];
cam_pts =[
    pts(1,:); pts(2,:); pts(3,:)
    pts(1,:); pts(3,:); pts(4,:)
    pts(1,:); pts(4,:); pts(5,:)
    pts(1,:); pts(5,:); pts(2,:)
   ]';
% cam_pts=[
%     -1 1 sz
%     0 0 0
%     1 -1 sz
%     1 1 sz
%     0 0 0
%     -1 -1 sz
%     1 -1 sz
%     1 1 sz
%     -1 1 sz
%     -1 -1 sz]';
cam_pts(4,:)=1;
c=[];
for i=1:m
    if ~isempty(opts.imgs)
        siz=size(opts.imgs{i});
    else
        siz=[1 1];
    end
    mx=max(siz);
    %scale matrix
    S=.02*[siz(2)/mx 0 0; 0 siz(1)/mx 0; 0 0 1];
    R=P(:,1:3,i);
    %fprintf('%g',det(R));
    t=P(:,4,i);
    %invers transformation from cam_i to the [I 0] coordinate
    g=[S*R' -R'*t];
    %g=[S*R t];
    c(:,end+1)=g*[0 0 0 1]';
    cam=double(g*cam_pts);
    if(isfield(opts,'cam_fill'))
        fill3(cam(1,:), cam(2,:), sz*cam(3,:), opts.cam_fill);
    end
    plot3(cam(1,:), cam(2,:), sz*cam(3,:), opts.cam_color);
    cnt=length(vertices);
    vertices=[vertices; [cam(1:3,:)' zeros(length(cam),3)]];
    edges=[edges; [cnt:cnt+length(cam)-2; cnt+1:cnt+length(cam)-1]' ];



    if ~isempty(opts.imgs)
        %--------- plot image on camera
        [X,Y,Z]=meshgrid(-1:2/(siz(2)-1):1, -1:2/(siz(1)-1):1, 1);
        X_=reshape(X,1,[]);
        Y_=reshape(Y,1,[]);
        Z_=reshape(Z,1,[]);
        A=[X_;Y_;Z_;ones(size(X_))];
        XX=g*A;
        X=reshape(XX(1,:),size(X));
        Y=reshape(XX(2,:),size(Y));
        Z=reshape(XX(3,:),size(Z));
        [img,map] = gray2ind( fliplr( imresize(opts.imgs{i},0.3) ) ,256);
        warp(X, Y, Z, (img));
        colormap(map);
        
        %----------------
        fvc=surf2patch(X, Y, Z, img);
        cnt=length(vertices);
        if size(img,3)==1
            color=double([fvc.facevertexcdata fvc.facevertexcdata fvc.facevertexcdata]);
        else
            color=double(fvc.facevertexcdata);
        end
%         vertices=[vertices; fvc.vertices color];
%         faces=[faces; fvc.faces+cnt-1];
    end
end
axis equal;
plot3(c(1,:), c(2,:), sz*c(3,:), ':k');
% text(c(1,:), c(2,:), sz*c(3,:), cellstr(num2str((1:size(c,2))')));


if ~isempty(opts.ply_file)
    f_id = fopen(opts.ply_file,'w');
    fprintf(f_id, ['ply\n',...
        'format ascii 1.0\n',...
        'element vertex %d\n',...
        'property float x\n',...
        'property float y\n',...
        'property float z\n',...
        'property uchar red\n',...
        'property uchar green\n',...
        'property uchar blue\n',...
        'element face %d\n',...
        'property list uchar int vertex_index\n',...
        'element edge %d\n',...
        'property int vertex1\n',...
        'property int vertex2\n',...
        'property uchar red\n',...
        'property uchar green\n',...
        'property uchar blue\n',...
        'end_header\n'],length(vertices),0,length(edges));%length(faces),length(edges));
    fprintf(f_id, '%f %f %f %g %g %g\n',vertices' );
    %fprintf(f_id, '4 %d %d %d %d\n',faces');
    fprintf(f_id, '%d %d %g %g %g\n',[edges ones(length(edges),3)]');
    fclose(f_id);
end
end
%
function opts=pars_options(varargin)
%--- defaults
opts.cam_color='r';
opts.K=eye(3);
opts.imgs=[];
opts.ply_file='';

%--- start parse
if nargin > 0
    i = 1;
    while i <= nargin
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'cam_color'
                    opts.cam_color=varargin{i+1};
                    i = i + 1;
                case 'cam_fill'
                    opts.cam_fill=varargin{i+1};
                    i = i + 1;
                case 'k'
                    opts.K=varargin{i+1};
                    i = i + 1;
                case 'imgs'
                    opts.imgs=varargin{i+1};
                    i = i + 1;
                case 'ply_file'
                    opts.ply_file=varargin{i+1};
                    i = i + 1;
            end
        end
        i = i + 1;
    end
end
end
