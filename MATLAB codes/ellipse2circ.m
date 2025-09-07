function [deltav_circ, r_circ, e_circ, i_circ, OM_circ, om_circ, theta_i_circ] = ellipse2circ(mu,a,e,i,OM,theta_circ)

% Function that circularize an elliptical orbit
%
% [deltav_circ, r_circ, e_circ, i_circ, OM_circ, om_circ, theta_i_circ] = ellipse2circ(mu, a, e, i, OM, theta_circ)
% 
% Input arguments:
% -------------------------------------------------------------------------
% mu           [1x1]     gravitational parameter                 [km^3/s^2]
% a            [1x1]     semi-major axis                         [km]
% e            [1x1]     eccentricity                            [-]
% i            [1x1]     inclination                             [rad]
% OM           [1x1]     RAAN                                    [rad]
% theta_circ   [1x1]     true anomaly of the manoeuvre           [rad]
%
% Output arguments:
% -------------------------------------------------------------------------
% deltav_circ  [1x1]     velocity difference in manoeuvre point  [km/s]
% r_circ       [1x1]     radius of circular orbit                [km]
% e_circ       [1x1]     eccentricity of circular orbit          [-]
% i_circ       [1x1]     inclination of circular orbit           [rad]
% OM_circ      [1x1]     RAAN of circular orbit                  [rad]
% om_circ      [1x1]     pericenter anomaly of circular orbit    [rad]
% theta_i_circ [1x1]     initial true anomaly of circular orbit  [rad]
% -------------------------------------------------------------------------

% calculate radius of the manoeuvre point:
r_circ = a*(1-e^2)/(1+e*cos(theta_circ));

% velocity on elliptical orbit at given manoeuvre true anomaly:
v_transv_ellipse = sqrt(mu/(a*(1-e^2))) * (1+e*cos(theta_circ));
v_radial_ellipse = sqrt(mu/(a*(1-e^2))) * e * sin(theta_circ);

% velocity on circular orbit at given manoeuvre true anomaly:
v_transv_circ = sqrt(mu/r_circ);
v_radial_circ = 0;

% velocity difference to circularize elliptical orbit:
deltav_transv = v_transv_circ - v_transv_ellipse;
deltav_radial = v_radial_circ - v_radial_ellipse;
deltav_circ = sqrt(deltav_transv^2 + deltav_radial^2);

% circular orbit parameters:
e_circ = 0;
i_circ = i;
OM_circ = OM;
om_circ = 0; % circular orbit has a 0 eccentricity
theta_i_circ = theta_circ;

end