function [deltav_cap, theta_pmv1, theta_pmv2] = changePericenterArg(mu,a_i,e_i,om_cp,om_f,theta)

% Function to change pericenter anomaly
%
% [deltav_cap, theta_pmv1, theta_pmv2] = changePericenterArg(mu, a_i, e_i, om_cp, om_f, theta)
%
% Input arguments:
% -------------------------------------------------------------------------
% mu          [1x1]   gravitational parameter                    [km^3/s^2]
% a_i         [1x1]   semi-major axis of the initial orbit       [km]
% e_i         [1x1]   eccentricity of the initial orbit          [-]
% om_cp       [1x1]   pericenter anomaly of the plane            [rad]
%                     changed orbit
% om_f        [1x1]   pericenter anomaly of the final orbit      [rad]
% theta       [1x1]   true anomaly of the starting point         [rad]
%                     in initial orbit
%
% Output arguments:
% -------------------------------------------------------------------------
% deltav_cap  [1x1]   velocity difference to change              [km/s]
%                     pericenter anomaly
% theta_pmv1  [1x1]   anomaly of the first possible manoeuvre    [rad]
%                     point to change pericenter anomaly
% theta_pmv2  [1x1]   anomaly of the second possible manoeuvre   [rad]
%                     point to change pericenter anomaly
% -------------------------------------------------------------------------

% calculate pericenter anomaly difference:
% pericenter anomaly of plane changed orbit is higher than the final one:
if om_cp > om_f
    dom = (2*pi+om_f) - om_cp;
    % comparing theta with dom to determine the anomaly of the two
    % manoeuvre points:
    if theta <=(2*pi+dom/2) && theta > (dom/2 + pi)
        theta_pmv1 = dom/2;
        theta_pmv2 = dom/2 + pi;
    elseif theta >= dom/2 && theta < (dom/2 + pi)
        theta_pmv1 = dom/2 + pi;
        theta_pmv2 = dom/2;
    end
    % calculate the velocity difference to change pericenter anomaly:
    deltav_cap = 2*abs(sqrt(mu/(a_i*(1-e_i^2)))*e_i*sin(dom/2));

% pericenter anomaly of plane changed orbit is lower than the final one:
elseif om_cp < om_f
    dom = om_f-om_cp;
    % comparing theta with dom to determine the anomaly of the two
    % manoeuvre points:
    if theta <= ((dom)/2 + pi) && theta > (dom)/2
        theta_pmv1 = (dom)/2 + pi;
        theta_pmv2 = (dom)/2;
    end
    if theta > ((dom)/2 + pi) && theta <= (2*pi+(dom)/2)
        theta_pmv1 = (dom)/2;
        theta_pmv2 = (dom)/2 + pi;
    end
    % calculate the velocity difference to change pericenter anomaly:
    deltav_cap = 2*sqrt(mu/(a_i*(1-e_i^2)))*e_i*sin((dom)/2);
end