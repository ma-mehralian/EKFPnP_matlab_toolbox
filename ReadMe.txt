%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this toolbox is an addition to the toolbox provided by the authors of
% MLPnP, CEPPnP and OPnP
% we extended it to show the use of EKFPnP
%
% Copyright (C) <2018>  <MA.Mehralian>
%
%     email: ma.mehralian@gmail.com
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Scripts to produce the synthetic results:
% main_ordinary_3d.m        
% main_planar.m          
% main_ordinary_3d_sigma.m   
% main_planar_sigma.m        
% main_ordinary_3d_time.m    
% main_planar_time.m	   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Acknowledgements
% This toolbox is based on the toolbox provided in:

% Urban, S.; Leitloff, J.; Hinz, S..
% MLPNP - A REAL-TIME MAXIMUM LIKELIHOOD SOLUTION TO THE PERSPECTIVE-N-POINT PROBLEM.
% ISPRS Annals of Photogrammetry, Remote Sensing & Spatial Information Sciences, 2016
% https://github.com/urbste/MLPnP_matlab

% Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer.
% Leveraging Feature Uncertainty in the PnP Problem. 
% In Proceedings of BMVC, 2014. 

% Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer.
% Very Fast Solution to the PnP Problem with Algebraic Outlier Rejection. 
% In Proceedings of CVPR, 2014.
%
% Y. Zheng, Y. Kuang, S. Sugimoto, K. Astro ?m, and M. Oku-
% tomi. Revisiting the pnp problem: A fast, general and opti-
% mal solution. In ICCV, pages 4321?4328, 2013.
%
% Y. Zheng, S. Sugimoto, and M. Okutomi. Aspnp: An accu- 
% rate and scalable solution to the perspective-n-point problem. 
% Trans. on Information and Systems, 96(7):1525?1535, 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



