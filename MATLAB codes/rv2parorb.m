function [a,e,i,OM,om,theta] = rv2parorb (rr, vv, mu)

% Trasformation from Cartesian state to orbital elements
%
% [a, e, i, OM, om, theta] = rv2parorb(r, v, mu)
%
% Input arguments:
% -------------------------------------------------------------------------
% rr        [3x1]   position vector                             [km]
% vv        [3x1]   velocity vector                             [km/s]
% mu        [1x1]   gravitational parameter                     [km^3/s^2]
%
% Output arguments:
% -------------------------------------------------------------------------
% a         [1x1]   semi-major axis                             [km]
% e         [1x1]   eccentricity                                [-]
% i         [1x1]   inclination                                 [rad]
% OM        [1x1]   RAAN                                        [rad]
% om        [1x1]   pericenter anomaly                          [rad]
% theta     [1x1]   true anomaly                                [rad]
% -------------------------------------------------------------------------


% step1: declare r and v, norms of the vectors rr and vv
r = norm(rr);
v = norm(vv);
% declare cartesian versors respectively of x,y,z axis
ivers = [1 0 0]';
jvers = [0 1 0]';
kvers = [0 0 1]';

% step2: calculate semi-major axis through specific orbital energy
% equation:
a = 1/(2/r - v.^2 / mu);

% step3: calculate specific angular momentum vector through vectorial
% products between position and velocity vectors, and its norm:
hvet = cross(rr,vv);
h = norm(hvet);

% step4: calculate eccentricity vector and its norm:
evet = (cross(vv,hvet)/mu) - rr/r;
e = norm(evet);

% step5: calculate orbit inclination:
i = acos((hvet'*kvers)./h);

% step6: calculate orbital node versor N:
if i == 0
    disp('equatorial orbit');
    OM = 0;
    Nvers = ivers;
else
    Nvers = cross(kvers,hvet)/norm(cross(kvers,hvet));

    % step7: calculate RAAN:
    if Nvers' * jvers >= 0
        OM = acos(Nvers' * ivers);
    else
        OM = 2*pi - acos(Nvers' * ivers);
    end
end

% step8: calculate pericenter anomaly:
if evet' * kvers >= 0
    om = acos((Nvers' * evet)/e);
else
    om = 2*pi - acos((Nvers' * evet)/e);
end

% step9: calculate true anomaly:
if vv' * rr >= 0
    theta = acos((rr' * evet)/(r*e));
else
    theta = 2*pi - acos((rr'*evet)/(r*e));
end

end