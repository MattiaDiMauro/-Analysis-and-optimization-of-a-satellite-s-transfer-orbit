function [deltav_1, deltav_2, delta_t_bitan, a_bitan, e_bitan] = bitangentTransfer(mu,a_i,e_i,a_f,e_f,type)

% Function to apply a bitangent transfer from an apsidal point of the first
% orbit to an apsidal point of the second orbit
%
% [deltav_1, deltav_2, delta_t_bitan, a_bitan, e_bitan] = bitangentTransfer(mu, a_i, e_i, a_f, e_f, type)
%
% Input arguments:
% -------------------------------------------------------------------------
% mu          [1x1]   gravitational parameter                    [km^3/s^2]
% a_i         [1x1]   semi-major axis of the initial orbit       [km]
% e_i         [1x1]   eccentricity of the initial orbit          [-]
% a_f         [1x1]   semi-major axis of the final orbit         [km]
% e_f         [1x1]   eccentricity of the final orbit            [-]
% type        [1x1]   type of bitangent transfer                 [-]
% type = 1 means a bitangent transfer that starts from pericenter of the
% first orbit and ends in apocenter of the second orbit
% type = 2 means a bitangent transfer that starts from apocenter of the
% first orbit and ends in pericenter of the second orbit
% 
% Output arguments:
% -------------------------------------------------------------------------
% deltav_1    [1x1]   velocity difference in the first           [km/s]
%                     apsidal point
% deltav_2    [1x1]   velocity difference in the second          [km/s]
%                     apsidal point
% delta_t_bitan [1x1]  time lapse spent on bitangent orbit       [s]
% a_bitan      [1x1]  semi-major axis of bitangent orbit         [km]
% e_bitan      [1x1]  eccentricity of the bitangent orbit        [km]
% -------------------------------------------------------------------------

% calculate semi-latus rectum of first and second orbit:
p_i = a_i * (1-e_i^2);
p_f = a_f * (1-e_f^2);

% calculate pericenter and apocenter radius of the bitangent orbit:
if type == 1
    r_p_bitan = p_i/(1+e_i);
    r_a_bitan = p_f/(1-e_f);
end

if type == 2
    r_a_bitan = p_i/(1-e_i);
    r_p_bitan = p_f/(1+e_f);
end

% calculate semi-major axis, eccentricity and semi-latus rectum of the
% bitangent orbit:
a_bitan = (r_p_bitan + r_a_bitan)/2;
e_bitan = (r_a_bitan - r_p_bitan)/(r_a_bitan + r_p_bitan);
p_bitan = a_bitan*(1-e_bitan^2);

if type == 1
    % calculate apsidal velocities:
    v_p_i = sqrt(mu/p_i)*(1+e_i);
    v_p_bitan = sqrt(mu/p_bitan)*(1+e_bitan);
    
    v_a_bitan = sqrt(mu/p_bitan)*(1-e_bitan);
    v_a_f = sqrt(mu/p_f)*(1-e_f);
    
    % calculate velocity difference of the two apsidal manoeuvre points:
    deltav_1 = v_p_bitan - v_p_i;
    deltav_2 = v_a_f - v_a_bitan;
end

if type == 2
    % calculate apsidal velocities:
    v_a_i = sqrt(mu/p_i)*(1-e_i);
    v_a_bitan = sqrt(mu/p_bitan)*(1-e_bitan);

    v_p_bitan = sqrt(mu/p_bitan)*(1+e_bitan);
    v_p_f = sqrt(mu/p_f)*(1+e_f);

    % calculate velocity difference of the two apsidal manoeuvre points:
    deltav_1 = v_a_bitan - v_a_i;
    deltav_2 = v_p_f - v_p_bitan;
end

% calculate time lapse spent on bitangent orbit:
delta_t_bitan = pi * sqrt(a_bitan^3 / mu);

end