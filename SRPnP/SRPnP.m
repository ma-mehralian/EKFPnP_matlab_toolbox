% Author :  Ping Wang                                                        
% Contact:  pingwangsky@gmail.com 
% This programe is implemented in matlab 2018a
% License:  Copyright (c) 2019 Ping Wang, All rights reserved       
% Address:  College of Electrical and Information Engineering, Lanzhou University of Technology              
% My site:  https://sites.google.com/view/ping-wang-homepage  

function [R t] = SRPnP( XX,xx )
%Thie version return all mimima to the end solution;
n= length(xx);
XXw= XX;

A=zeros(2*n,3);
B1=zeros(2*n,10);
uuj=xx(1,:)'; vvj=xx(2,:)';
Wx=XXw(1,:)';Wy=XXw(2,:)';Wz=XXw(3,:)';
A(1:2:end,:)=[ones(n,1),zeros(n,1),-uuj];
A(2:2:end,:)=[zeros(n,1),ones(n,1),-vvj];
B1(1:2:end,:)=[-Wx+Wz.*uuj,2*Wy.*uuj,-2*Wz-2*Wx.*uuj,2*Wy,-Wx-Wz.*uuj,...
               -2*Wy,-2*Wz+2*Wx.*uuj,Wx-Wz.*uuj,2*Wy.*uuj,Wz.*uuj+Wx;];
B1(2:2:end,:)=[-Wy+Wz.*vvj,2*Wy.*vvj+2*Wz,-2*Wx.*vvj,-2*Wx,Wy-Wz.*vvj,...
               -2*Wx,2*Wx.*vvj,-Wy-Wz.*vvj,-2*Wz+2*Wy.*vvj,Wz.*vvj+Wy];
C1=((A.'*A)\A.')*B1; 
E1=A*C1-B1;
G1=E1.'*E1;

%normalized;
xxv= [xx; ones(1,n)];
for i=1:n
    xxv(:,i)= xxv(:,i)/norm(xxv(:,i));
end

% selecting an edge $P_{i1}P_{i2}$ whose projection length
% is the longest in image plane.
lmin= inf;
i1= 0; i2= 0;
for i= 1:n-1
    for j= i+1:n
        l= xxv(1,i)*xxv(1,j)+xxv(2,i)*xxv(2,j)+xxv(3,i)*xxv(3,j);
        if l < lmin
            i1= i;
            i2= j;
            lmin= l;
        end
    end
end

% calculating the rotation matrix of $O_aX_aY_aZ_a$.
p1= XX(:,i1);
p2= XX(:,i2);
p0= (p1+p2)/2;
x= p2-p0; x= x/norm(x);

if abs([0 1 0]*x) < abs([0 0 1]*x)
    z= xcross(x,[0; 1; 0]); z= z/norm(z);
    y= xcross(z, x); y= y/norm(y);
else
    y= xcross([0; 0; 1], x); y= y/norm(y);
    z= xcross(x,y); z= z/norm(z);
end
Ro= [x y z].';

% transforming the reference points form orignial object space 
% to the new coordinate frame  $O_aX_aY_aZ_a$.
XX= Ro*(XX-repmat(p0,1,n));

% Dividing the n-point set into (n-2) 3-point subsets
% and setting up the P3P equations

v1= xxv(:,i1);
v2= xxv(:,i2);

cg1= v1.'*v2;
sg1= sqrt(1-cg1^2);
D1= norm(XX(:,i1)-XX(:,i2));
D4= zeros(n-2,5);

j= 0;
for i= 1:n
    if i == i1 || i == i2
        continue;
    end
    j= j+1;

    vi= xxv(:,i);
    cg2= v1.'*vi;
    cg3= v2.'*vi;
    sg2= sqrt(1-cg2^2);
    D2= norm(XX(:,i1)-XX(:,i));
    D3= norm(XX(:,i)-XX(:,i2));
        
    % get the coefficients of the P3P equation from each subset.
     D4(j,:)= getp3p(cg1,cg2,cg3,sg1,sg2,D1,D2,D3);
     
end

% get the 7th order polynomial, the deviation of the cost function.
D7= zeros(1,8);
for i= 1:n-2
   D7= D7+ getpoly7(D4(i,:));
end

% retriving the local minima of the cost function.
t2s= roots(D7);

maxreal= max(abs(real(t2s)));
t2s(abs(imag(t2s))/maxreal > 0.001)= [];
t2s= real(t2s);

D6= (7:-1:1).*D7(1:7);
F6= polyval(D6,t2s);
t2s(F6 <= 0)= [];

if isempty(t2s)
    fprintf('no solution!\n');
    return
end

% calculating the camera pose from each local minimum.
m= length(t2s);
Px=XX(1,:)';Py=XX(2,:)';Pz=XX(3,:)';
count=1;
for i=1:m
    t2= t2s(i);
    % calculating the rotation matrix
    d2= cg1+t2;
    x= v2*d2- v1; x= x/norm(x);
    
    if abs([0 1 0]*x) < abs([0 0 1]*x)
        z= xcross(x,[0; 1; 0]); z= z/norm(z);
        y= xcross(z, x); y= y/norm(y);
    else
        y= xcross([0; 0; 1], x); y= y/norm(y);
        z= xcross(x,y); z= z/norm(z);
    end
    Rx= [x y z];
    
    % calculating p, q;
    B=zeros(2*n,3);
    B(1:2:end,:)=[Rx(3,2)*Py.*uuj+Rx(3,3)*Pz.*uuj-Rx(1,2)*Py-Rx(1,3)*Pz,...
                  Rx(3,3)*Py.*uuj-Rx(3,2)*Pz.*uuj-Rx(1,3)*Py+Rx(1,2)*Pz,...
                  Rx(3,1)*Px.*uuj-Rx(1,1)*Px];
    B(2:2:end,:)=[Rx(3,2)*Py.*vvj+Rx(3,3)*Pz.*vvj-Rx(2,2)*Py-Rx(2,3)*Pz,...
                  Rx(3,3)*Py.*vvj-Rx(3,2)*Pz.*vvj-Rx(2,3)*Py+Rx(2,2)*Pz,...
                  Rx(3,1)*Px.*vvj-Rx(2,1)*Px]; 
    C=((A.'*A)\A.')*B;   
    
    E=A*C-B;  
    G=E.'*E;  
    
    % get g from G;
     g11=G(1,1);
     g12=G(1,2);
     g22=G(2,2);
     g23=G(2,3);
     g13=G(1,3);
     
     AA=G(2,2)-G(1,1);
     F4=4*g12^2+AA^2;
     F3=4*g12*g23-2*g13*AA;
     F2=g23^2-4*g12^2-AA^2+g13^2;
     F1=-2*g12*g23+2*g13*AA;
     F0=g12^2-g13^2;
    
     c= roots([F4,F3,F2,F1,F0]);

     maxreal1= max(abs(real(c)));
     c(abs(imag(c))/maxreal1 > 0.001)= [];

    c= real(c);
    c=c.';
    s=(2*g12*c.^2+g23*c-g12)./((g11-g22)*c+g13);
    n1=length(c);
    for j=1:n1   
        ss=[c(j);s(j);1];
        ts=C*ss;
        %calculating the camera pose by 3d alignment
        xi= XX(1,:); yi= XX(2,:); zi= XX(3,:);
        r1=Rx(1,1);r2=Rx(1,2);r3=Rx(1,3);
        r4=Rx(2,1);r5=Rx(2,2);r6=Rx(2,3);
        r7=Rx(3,1);r8=Rx(3,2);r9=Rx(3,3);
        XXcs=[r1*xi+(r2*c(j)+r3*s(j))*yi+(-r2*s(j)+r3*c(j))*zi+ts(1);
              r4*xi+(r5*c(j)+r6*s(j))*yi+(-r5*s(j)+r6*c(j))*zi+ts(2);
              r7*xi+(r8*c(j)+r9*s(j))*yi+(-r8*s(j)+r9*c(j))*zi+ts(3)];
        XXc= zeros(size(XXcs));  
        for j1= 1:n
            XXc(:,j1)= xxv(:,j1)*norm(XXcs(:,j1));
        end      
        [R1,~]= calcampose(XXc,XXw);
        
        % refine solution by using single Gauss-Newton method; 
        solution = RefineGaussNewton(Cayley(R1),E1,G1);
        
        s1=solution(1);
        s2=solution(2);
        s3=solution(3);
        Rw=[1+s1^2-s2^2-s3^2,2*s1*s2-2*s3,2*s2+2*s1*s3;
           2*s3+2*s1*s2,1-s1^2+s2^2-s3^2,2*s2*s3-2*s1;
           2*s1*s3-2*s2,2*s1+2*s2*s3,1-s1^2-s2^2+s3^2];
        w=[1,s1,s2,s3,s1^2,s1*s2,s1*s3,s2^2,s2*s3,s3^2].';

        factor=1/(1+s1^2+s2^2+s3^2);
        R2=Rw*factor;
        t2=C1*w*factor;
              
        R(:,:,count)=R2;
        t(:,count)=t2;
        
        count=count+1;
    end         
end

end



