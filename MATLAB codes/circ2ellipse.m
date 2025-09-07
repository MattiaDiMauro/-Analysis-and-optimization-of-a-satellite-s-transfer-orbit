function [deltav_1] = circ2ellipse(mu, r_circ, a, e, maintain)

% Function to change shape of a circular orbit into an elliptical orbit
%
% [deltav_1] = circ2ellipse(mu, r_circ, a, e, maintain)
%
% Input arguments:
% -------------------------------------------------------------------------
% mu           [1x1]     gravitational parameter                 [km^3/s^2]
% r_circ       [1x1]     circular orbit radius                   [km]
% a            [1x1]     semi-major axis of final orbit          [km]
% e            [1x1]     eccentricity of final orbit             [-]
% mantain      [1x1]     variable to understand if               [-]
%                        circular orbit radius will
%                        become periapsis or apoapsis
%
% Output arguments:
% -------------------------------------------------------------------------
% deltav_1     [1x1]     velocity difference to                  [km/s]
%                        change shape of the orbit
% -------------------------------------------------------------------------

% maintain = 1 means apoapsis will have the same value of circular orbit
% radius
if maintain == 1
    theta_mvr = pi;
% maintain = 2 means periapsis will have the same value of circular orbit
% radius
elseif maintain == 2
    theta_mvr = 0;
end

% velocity on circular orbit at given manoeuvre true anomaly:
v_transv_circ = sqrt(mu/r_circ);

% velocity on elliptical orbit at given manoeuvre true anomaly:
v_transv_ellipse = sqrt(mu/(a*(1-e^2))) * (1+e*cos(theta_mvr));

% velocity difference in first manoeuvre point:
deltav_transv = v_transv_ellipse - v_transv_circ;
deltav_1 = deltav_transv;

end