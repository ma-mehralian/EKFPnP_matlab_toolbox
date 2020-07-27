% Author :  Ping Wang                                                        
% Contact:  pingwangsky@gmail.com 
% This programe is implemented in matlab 2018a
% License:  Copyright (c) 2019 Ping Wang, All rights reserved       
% Address:  College of Electrical and Information Engineering, Lanzhou University of Technology              
% My site:  https://sites.google.com/view/ping-wang-homepage  

function solution= RefineGaussNewton(solution,E,G1)
    %refine the solution by using one step Gauss Newton method
    s1=solution(1); s2=solution(2); s3=solution(3);
    w=[1,s1,s2,s3,s1^2,s1*s2,s1*s3,s2^2,s2*s3,s3^2].';
    obj_pre = w.'*G1*w;
    s1=solution(1); s2=solution(2); s3=solution(3);
    w=[1,s1,s2,s3,s1^2,s1*s2,s1*s3,s2^2,s2*s3,s3^2].';
    %Jacobian matrix
    Jac=[ 0,     0,     0;
          1,     0,     0;
          0,     1,     0;
          0,     0,     1;
       2*s1,     0,     0;
         s2,    s1,     0;
         s3,     0,    s1;
          0,  2*s2,     0;
          0,    s3,    s2;
          0,     0,   2*s3 ];
     Fk=E*w;  
     Jk=E*Jac;
    
    solution_temp = solution;
    %increment
    dk=-(Jk.'*Jk)\(Jk.'*Fk);
    %update parameter
    solution_temp = solution_temp + dk; 
    s1=solution_temp(1); s2=solution_temp(2); s3=solution_temp(3);
    w=[1,s1,s2,s3,s1^2,s1*s2,s1*s3,s2^2,s2*s3,s3^2].';
    %evaluate the error of objection; 
    obj_cur = w.'*G1*w;
    if obj_cur<obj_pre 
        solution=solution_temp;
    end
end

