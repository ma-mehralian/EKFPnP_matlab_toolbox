% Author :  Ping Wang                                                        
% Contact:  pingwangsky@gmail.com 
% This programe is implemented in matlab 2018a
% License:  Copyright (c) 2019 Ping Wang, All rights reserved       
% Address:  College of Electrical and Information Engineering, Lanzhou University of Technology              
% My site:  https://sites.google.com/view/ping-wang-homepage    

function F7= getpoly7(F)

F7= [4*F(1)^2;
7*F(2)*F(1);
6*F(3)*F(1)+3*F(2)^2;
5*F(4)*F(1)+5*F(3)*F(2);
4*F(5)*F(1)+4*F(4)*F(2)+2*F(3)^2;
3*F(5)*F(2)+3*F(4)*F(3);
2*F(5)*F(3)+F(4)^2;
F(5)*F(4)].';

end