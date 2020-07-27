% Author :  Ping Wang                                                        
% Contact:  pingwangsky@gmail.com 
% This programe is implemented in matlab 2018a
% License:  Copyright (c) 2019 Ping Wang, All rights reserved       
% Address:  College of Electrical and Information Engineering, Lanzhou University of Technology              
% My site:  https://sites.google.com/view/ping-wang-homepage  

function s= Cayley(R)

A= R.';

q4= sqrt(1+A(1,1)+A(2,2)+A(3,3))/2;

if q4 > 0.01
	q1= (A(3,2)-A(2,3))/q4/4;
	q2= (A(1,3)-A(3,1))/q4/4;
	q3= (A(2,1)-A(1,2))/q4/4;
else
	q1= sqrt(1+A(1,1)-A(2,2)-A(3,3))/2;
	q2= (A(1,2)+A(2,1))/q1/4;
	q3= (A(1,3)+A(3,1))/q1/4;
	q4= (A(3,2)-A(2,3))/q1/4;
end

s= -[q1; q2; q3]/q4;

return