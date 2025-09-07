function [deltav_cp,om_cp,theta_cp1,theta_cp2] = changeOrbitalPlane(mu,a_i,e_i,i_i,OM_i,om_i,i_f,OM_f,theta)

% Function to change orbital plane
%
% [deltav_cp, om_cp, theta_cp1, theta_cp2] = changeOrbitalPlane(mu,a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, theta)
%
% Input arguments:
% -------------------------------------------------------------------------
% mu          [1x1]   gravitational parameter                    [km^3/s^2]
% a_i         [1x1]   semi-major axis of the initial orbit       [km]
% e_i         [1x1]   eccentricity of the initial orbit          [-]
% i_i         [1x1]   inclination of the initial orbit           [rad]
% OM_i        [1x1]   RAAN of the initial orbit                  [rad]
% om_i        [1x1]   pericenter anomaly of the initial orbit    [rad]
% i_f         [1x1]   inclination of the final orbit             [rad]
% OM_f        [1x1]   RAAN of the final orbit                    [rad]
% theta       [1x1]   true anomaly of the starting point         [rad]
%                     in initial orbit
%
% Output arguments:
% -------------------------------------------------------------------------
% deltav_cp   [1x1]   velocity difference to change plane        [km/s]
% om_cp       [1x1]   pericenter anomaly of the plane            [rad]
%                     changed orbit
% theta_cp1   [1x1]   anomaly of the first possible              [rad]
%                     manoeuvre point to change plane
% theta_cp2   [1x1]   anomaly of the second possible             [rad]
%                     manoeuvre point to change plane
% -------------------------------------------------------------------------

% solve spherical triangle:
cos_alpha = -cos(pi-i_f)*cos(i_i)+sin(pi-i_f)*sin(i_i)*cos(OM_f-OM_i);
alpha = acos(cos_alpha);

sin_u2=(sin(i_i)*sin(OM_f-OM_i)/sin(alpha));
sin_u1=(sin(pi-i_f)*sin(OM_f-OM_i)/sin(alpha));

cos_u2=(cos(i_i)+cos(alpha)*cos(pi-i_f))/(sin(alpha)*sin(pi-i_f));
cos_u1=(cos(pi-i_f)+cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));

u2=atan2(sin_u2,cos_u2);
u1=atan2(sin_u1,cos_u1);

theta_cp=u1-om_i;

% calculate pericenter anomaly of the plane changed orbit:
om_cp=u2-theta_cp;
if om_cp < 0
    om_cp = 2*pi-abs(om_cp);
end

% comparing theta_cp obtained from the spherical triangle with theta 
% to determine which is the first possible anomaly to change plane:
if theta > theta_cp
    theta_cp1 = theta_cp+pi;
    theta_cp2 = theta_cp;
elseif theta > (theta_cp+pi)
    theta_cp1 = theta_cp;
    theta_cp2 = theta_cp+pi;
end

% calculate semi-latus rectum of the orbit:
p=a_i*(1-e_i^2);

% calculate transverse velocity in point whose anomaly is theta_cp1
v_theta=sqrt(mu/p)*(1+e_i*cos(theta_cp1));

% calculate velocity difference to change plane:
deltav_cp=2*v_theta*sin(alpha/2);

end