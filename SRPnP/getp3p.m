% Author :  Ping Wang                                                        
% Contact:  pingwangsky@gmail.com 
% This programe is implemented in matlab 2018a
% License:  Copyright (c) 2019 Ping Wang, All rights reserved       
% Address:  College of Electrical and Information Engineering, Lanzhou University of Technology              
% My site:  https://sites.google.com/view/ping-wang-homepage  

function B = getp3p(l1,l2,A5,C1,C2,D1,D2,D3)
A1= (D2/D1)^2;
A2= A1*C1^2-C2^2;
A3= l2*A5-l1;
A4= l1*A5-l2;
A6= (D3^2-D1^2-D2^2)/(2*D1^2);
A7= 1-l1^2-l2^2+l1*l2*A5+A6*C1^2;

B= [A6^2-A1*A5^2, 2*(A3*A6-A1*A4*A5), A3^2+2*A6*A7-A1*A4^2-A2*A5^2,...
    2*(A3*A7-A2*A4*A5), A7^2-A2*A4^2];

end

