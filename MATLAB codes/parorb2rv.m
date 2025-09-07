function [rr,vv] = parorb2rv( a,e,i,OM,om,theta,mu)

% Trasformation from orbital elements to Cartesian coordinates
%
% [rr,vv] = parorb2rv(a, e, i, OM, om, theta, mu)
%
% Input arguments:
% -------------------------------------------------------------------------
% a         [1x1]   semi-major axis                             [km]
% e         [1x1]   eccentricity                                [-]
% i         [1x1]   inclination                                 [rad]
% OM        [1x1]   RAAN                                        [rad]
% om        [1x1]   pericenter anomaly                          [rad]
% theta     [1x1]   true anomaly                                [rad]
% mu        [1x1]   gravitational parameter                     [km^3/s^2]
%
% Output arguments:
% -------------------------------------------------------------------------
% rr        [3x1]   position vector                             [km]
% vv        [3x1]   velocity vector                             [km/s]
% -------------------------------------------------------------------------

% step1: calculate semi-latus rectum:
p = a*(1-e^2);

% step2: calculate position vector expressed in perifocal system:
r_vet_PF = [p/(1+e*cos(theta)) * cos(theta),...
    p/(1+e*cos(theta)) * sin(theta), 0]';

% step3: calculate radial and transverse components of velocity vector:
v_vet_r0k = [sqrt(mu/p) * e * sin(theta), ...
    sqrt(mu/p) *(1+e*cos(theta)), 0]';

R_theta = [cos(theta), sin(theta), 0; ...
    -sin(theta), cos(theta), 0; 0 0 1];

% step4: obtain velocity vector expressed in perifocal system through a
% rotation with matric R_theta:
v_vet_PF = R_theta' * v_vet_r0k;

% step5: declare rotation matrixes:
R_OM = [cos(OM) sin(OM) 0; ...
    -sin(OM) cos(OM) 0; 0 0 1];

R_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];

R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

% step6: global rotation matrix from perifocal system to equatorial
% geocentric system:
R_PF_to_GE = (R_om * R_i * R_OM)';

% step7: apply the rotation to perifocal system expressed vectors and
% obtain new ones in equatorial geocentrical system:
rr = R_PF_to_GE * r_vet_PF;
vv = R_PF_to_GE * v_vet_PF;

end