function [delta_t] = TOF(mu, a, e, theta_i, theta_f)

% Time of Flight
%
% [delta_t] = TOF(mu, a, e, theta_i, theta_f)
%
% Input arguments:
% -------------------------------------------------------------------------
% mu        [1x1]   gravitational parameter                    [km^3/s^2]
% a         [1x1]   semi-major axis                            [km]
% e         [1x1]   eccentricity                               [-]
% theta_i   [1x1]   initial true anomaly                       [rad]
% theta_f   [1x1]   final true anomaly                         [rad]
%
% Output arguments:
% -------------------------------------------------------------------------
% delta_t   [1x1]   time of flight                             [s]
%
% -------------------------------------------------------------------------


delta_t = 0;

% this is the case of a theta_f < theta_i:
if theta_f > 2*pi
    theta_f = theta_f - 2*pi;
    % we already know that we will make a full lap around orbit:
    delta_t = 2*pi*sqrt(a^3/mu);
end

% calculate initial eccentric anomaly:
E_i = 2*atan( tan(theta_i/2)/sqrt((1+e)/(1-e)));
if E_i < 0
    E_i = 2*pi - abs(E_i);
end

% calculate final eccentric anomaly:
E_f = 2*atan( tan(theta_f/2)/sqrt((1+e)/(1-e)));
if E_f < 0
   E_f = 2*pi - abs(E_f);
end

% calculate time to reach each initial and final true anomaly:
t_i = sqrt(a^3 / mu) * (E_i - e*sin(E_i));
t_f = sqrt(a^3 / mu) * (E_f - e*sin(E_f));

% time of flight to reach final true anomaly starting from initial true
% anomaly:
delta_t = delta_t + t_f - t_i;

end