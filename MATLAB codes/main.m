%% STANDARD TRANSFER (CP + CAP + BITAN a_i-p_f)
clear;
close all;
clc;

% gravitational parameter:
mu = 398600; % [km^3/s^2]

% initial orbit data:

% Cartesian coordinates[km]
x = -5829.6395;
y = -6946.1419;
z = 4067.7268;
rvet = [x,y,z]';

% velocity vector components[km/s]
vx = 1.9470;
vy = -4.5900;
vz = -3.8490;
vvet = [vx,vy,vz]';

% Keplerian orbital parameters of initial orbit:
[a_i,e_i,i_i,OM_i,om_i,theta_i] = rv2parorb(rvet, vvet, mu);


% final orbit data:

% Keplerian orbital parameters of final orbit:
a_f = 13790.0000; %[km]
e_f = 0.2719; %[-]
i_f = 1.0470; % [rad]
OM_f = 2.5590; % [rad]
om_f = 0.4614; % [rad]
theta_f = 0.2578; % [rad]

% Plane change manoeuvre:
[deltav_cp,om_cp,theta_cp1,theta_cp2] = changeOrbitalPlane(mu,a_i,e_i,i_i,OM_i,om_i,i_f,OM_f,theta_i);

% Pericenter anomaly adjustment:
[deltav_cap, theta_pmv1, theta_pmv2] = changePericenterArg(mu,a_i,e_i,om_cp,om_f,theta_cp1);

% Bitangent transfer:
type = 2;
[deltav_3, deltav_4, deltat_4, a_bitan, e_bitan] = bitangentTransfer(mu,a_i,e_i,a_f,e_f,type);

% total velocity difference:
deltav_tot = abs(deltav_cp)+abs(deltav_cap)+abs(deltav_3)+abs(deltav_4);

% Graphical representation:
% declare initial,final true anomaly and discretisation step for initial
% orbit:
th0_i = theta_i;
thf_i = theta_cp1;
dth_i = 0.01;
% initial orbit plot:
color = 1;
plotOrbit(a_i, e_i, i_i, OM_i, om_i, th0_i, thf_i, dth_i, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for plane
% changed orbit:
th0_cp = theta_cp1;
thf_cp = 2*pi+theta_pmv1;
dth_cp = 0.01;
% plane changed orbit plot:
color = 2;
plotOrbit(a_i, e_i, i_f, OM_f, om_cp, th0_cp, thf_cp, dth_cp, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for pericenter
% anomaly adjusted orbit:
th0_cap = 2*pi-theta_pmv1;
thf_cap = 3*pi;
dth_cap = 0.01;
% pericenter anomaly adjusted orbit plot:
color = 3;
plotOrbit(a_i, e_i, i_f, OM_f, om_f, th0_cap, thf_cap, dth_cap, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for bitangent
% transfer orbit:
th0_bitan = pi;
thf_bitan = 2*pi;
dth_bitan = 0.01;
% bitangent transfer orbit plot:
color = 4;
plotOrbit(a_bitan,e_bitan,i_f,OM_f,om_f,th0_bitan,thf_bitan, dth_bitan, mu,color);

% declare initial,final true anomaly and discretisation step for final
% orbit:
th0_f = 0;
thf_f = theta_f;
dth_f = 0.01;
% final orbit plot:
color = 5;
plotOrbit(a_f, e_f, i_f, OM_f, om_f, th0_f, thf_f, dth_f, mu,color);
hold on

% earth plot:
earth3D

legend('initial orbit','plane change point','plane changed orbit','pericenter anomaly adjustment point','pericenter anomaly adjusted orbit','initial orbit apocenter','bitangent orbit','final orbit pericenter','final orbit','target point')

% Calculate time lapses:
% time lapse from relapse point to plane change point:
[delta_t1] = TOF(mu,a_i,e_i,th0_i, thf_i);

% time lapse from plane change point to pericenter anomaly adjustment
% point:
[delta_t2] = TOF(mu,a_i,e_i,th0_cp, thf_cp);

% time lapse from pericenter anomaly adjustment point to bitangent transfer
% orbit:
[delta_t3] = TOF(mu,a_i,e_i,th0_cap, thf_cap);

% time spent on bitangent transfer orbit:
[delta_t4] = TOF(mu,a_bitan,e_bitan,th0_bitan, thf_bitan);

% time lapse from final orbit pericenter to final target point:
[delta_t5] = TOF(mu,a_f,e_f,th0_f, thf_f);

% calculate total time lapse:
delta_t_totale = delta_t1+delta_t2+delta_t3+delta_t4+delta_t5;

%% ALTERNATIVE TRANSFER1 ~ Step changed transfer (BITAN p_i-a_f + CP + CAP)
clear;
clc;
close all;

% gravitational parameter:
mu = 398600; % [km^3/s^2]

% initial orbit data:
% Cartesian coordinates[km]
x = -5829.6395;
y = -6946.1419;
z = 4067.7268;
rvet = [x,y,z]';

% velocity vector components[km/s]
vx = 1.9470;
vy = -4.5900;
vz = -3.8490;
vvet = [vx,vy,vz]';

% Keplerian orbital parameters of initial orbit:
[a_i,e_i,i_i,OM_i,om_i,theta_i] = rv2parorb(rvet, vvet, mu);


% final orbit data:

% Keplerian orbital parameters of final orbit:
a_f = 13790.0000; %[km]
e_f = 0.2719; %[-]
i_f = 1.0470; % [rad]
OM_f = 2.5590; % [rad]
om_f = 0.4614; % [rad]
theta_f = 0.2578; % [rad]

% Bitangent transfer:
type = 1;
[deltav_1, deltav_2, deltat_1, a_bitan, e_bitan] = bitangentTransfer(mu,a_i,e_i,a_f,e_f,type);

% Plane change manoeuvre:
[deltav_cp,om_cp,theta_cp1,theta_cp2] = changeOrbitalPlane(mu,a_f,e_f,i_i,OM_i,om_i,i_f,OM_f,pi);

% Pericenter anomaly adjustment:
[deltav_cap, theta_pmv1, theta_pmv2] = changePericenterArg(mu,a_f,e_f,om_cp,om_f,theta_cp1);

% total velocity difference:
deltav_totale = abs(deltav_cp)+abs(deltav_cap)+abs(deltav_1)+abs(deltav_2);

% Graphical representation:
% declare initial,final true anomaly and discretisation step for initial
% orbit:
th0_i = theta_i;
thf_i = 2*pi;
dth_i = 0.01;
% initial orbit plot:
color = 1;
plotOrbit(a_i, e_i, i_i, OM_i, om_i, th0_i, thf_i, dth_i, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for bitangent
% transfer orbit:
th0_bitan = 0;
thf_bitan = pi;
dth_bitan = 0.01;
% bitangent transfer orbit plot:
color = 2;
plotOrbit(a_bitan,e_bitan,i_i,OM_i,om_i,th0_bitan,thf_bitan, dth_bitan, mu,color);

% declare initial, final true anomaly and discretisation step for orbit
% that has the same shape of the final one but in the initial orbit plane:
th0_fplanei = pi;
thf_fplanei = theta_cp1;
dth_fplanei = 0.01;
% plot of the orbit that has the same shape of the final one but in the
% initial orbit plane:
color = 3;
plotOrbit(a_f,e_f,i_i,OM_i,om_i,th0_fplanei,thf_fplanei, dth_fplanei, mu,color);

% declare initial,final true anomaly and discretisation step for plane
% changed orbit:
th0_cp = theta_cp1;
thf_cp = 2*pi+theta_pmv1;
dth_cp = 0.01;
% plane changed orbit plot:
color = 4;
plotOrbit(a_f, e_f, i_f, OM_f, om_cp, th0_cp, thf_cp, dth_cp, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for pericenter
% anomaly adjusted orbit:
th0_cap = 2*pi-theta_pmv1;
thf_cap = 2*pi+theta_f;
dth_cap = 0.01;
% pericenter anomaly adjusted orbit plot:
color = 5;
plotOrbit(a_f, e_f, i_f, OM_f, om_f, th0_cap, thf_cap, dth_cap, mu,color);
hold on

% earth plot:
earth3D

legend('initial orbit','initial orbit pericenter','bitangent orbit','bitangent orbit apocenter','final orbit initial plane','plane change point','plane changed orbit','pericenter anomaly adjustment point','final orbit','target point')

% Calculate time lapses:
% time laps from relapse point to pericenter of initial orbit:
[delta_t1] = TOF(mu,a_i,e_i,th0_i, thf_i);

% time spent on bitangent transfer orbit:
[delta_t2] = TOF(mu,a_bitan,e_bitan,th0_bitan, thf_bitan);

% time lapse from apocenter of bitangent transfer orbit to plane change
% point:
[delta_t3] = TOF(mu,a_f,e_f,th0_fplanei, thf_fplanei);

% time lapse from plane change point to pericenter anomaly adjustment
% point:
[delta_t4] = TOF(mu,a_f,e_f,th0_cp, thf_cp);

% time lapse from pericenter anomaly adjustment point to final target
% point:
[delta_t5] = TOF(mu,a_f,e_f,th0_cap, thf_cap);

% calculate total time lapse:
delta_t_totale = delta_t1+delta_t2+delta_t3+delta_t4+delta_t5;

%% ALTERNATIVE TRANSFER2 ~ Time reduction
clear;
close all;
clc;

% gravitational parameter:
mu = 398600; % [km^3/s^2]

% initial orbit data:
% Cartesian coordinates[km]
x = -5829.6395;
y = -6946.1419;
z = 4067.7268;
rvet = [x,y,z]';

% velocity vector components[km/s]
vx = 1.9470;
vy = -4.5900;
vz = -3.8490;
vvet = [vx,vy,vz]';

% Keplerian orbital parameters of initial orbit:
[a_i,e_i,i_i,OM_i,om_i,theta_i] = rv2parorb(rvet, vvet, mu);


% final orbit data:

% Keplerian orbital parameters of final orbit:
a_f = 13790.0000; %[km]
e_f = 0.2719; %[-]
i_f = 1.0470; % [rad]
OM_f = 2.5590; % [rad]
om_f = 0.4614; % [rad]
theta_f = 0.2578; % [rad]

% Final orbit periapsis is included between initial orbit periapsis and
% apoapsis, so we can operate a tangent manoeuvre in the point of initial
% orbit that has a radius equal to final orbit periapsis to circularize the
% initial orbit:
p_i = a_i*(1-e_i^2); % semi-latus rectum of initial orbit [km]
r_p_f = a_f*(1-e_f); % final orbit periapsis [km]
thet_tang = acos((1/e_i)*((p_i/r_p_f)-1)); % true anomaly of initial 
% orbit where we operate the tangent manoeuvre [rad]
[deltav_circ, r_circ, e_circ, i_circ, OM_circ, om_circ, theta_i_circ] = ellipse2circ(mu,a_i,e_i,i_i,OM_i,thet_tang);

% Plane change manoeuvre for circular orbit:
[deltav_cp,om_cp,theta_cp1,theta_cp2] = changeOrbitalPlane(mu,r_circ,e_circ,i_circ,OM_circ,om_circ,i_f,OM_f,theta_i_circ);

% Tangent manoeuvre to change shape of circular orbit into the final orbit:
% we will maintain periapsis as r_circ:
maintain = 2;
[deltav_3] = circ2ellipse(mu, r_circ, a_f, e_f, maintain);

% total velocity difference:
deltav_tot = abs(deltav_circ)+abs(deltav_cp)+abs(deltav_3);

% Graphical representation:
% declare initial,final true anomaly and discretisation step for initial
% orbit:
th0_i = theta_i;
thf_i = thet_tang;
dth_i = 0.01;
% initial orbit plot:
color=1;
plotOrbit(a_i, e_i, i_i, OM_i, om_i, th0_i, thf_i, dth_i, mu, color);
hold on

% declare initial, final true anomaly and discretisation step for circular 
% orbit:
th0_circ = thet_tang+om_i;
thf_circ = theta_cp1;
dth_circ = 0.01;
% circular orbit plot:
color=2;
plotOrbit(r_circ, e_circ, i_circ, OM_circ, om_circ, th0_circ, thf_circ, dth_circ, mu, color);
hold on

% declare initial, final true anomaly and discretisation step for plane
% changed circular orbit:
th0_cp = theta_cp1+om_cp - 2*pi;
thf_cp = 2*pi+om_f;
dth_cp = 0.01;
% plane changed circular orbit plot:
color=3;
plotOrbit(r_circ, e_circ, i_f, OM_f, om_circ, th0_cp, thf_cp, dth_cp, mu, color);
hold on

% declare initial, final true anomaly and discretisation step for final
% orbit:
th0_f = 0;
thf_f = theta_f;
dth_f = 0.01;
% final orbit plot:
color=4;
plotOrbit(a_f, e_f, i_f, OM_f, om_f, th0_f, thf_f, dth_f, mu, color);
hold on

% earth plot:
earth3D

legend('initial orbit','circularization point','circular orbit','plane change point','circular plane changed orbit','tangent manoeuvre point','final orbit','target point')

% calculate time lapses:

% time lapse from relapse point to tangent manoeuvre point:
[delta_t1] = TOF(mu,a_i,e_i,th0_i, thf_i);

% time lapse from tangent manoeuvre point to plane change point: 
[delta_t2] = TOF(mu,r_circ,e_circ,th0_circ, thf_circ);

% time lapse from plane change point to second tangent manoeuvre point:
[delta_t3] = TOF(mu,r_circ,e_circ,th0_circ, thf_circ);

% time lapse from second tangent manoeuvre point to final target point:
[delta_t4] = TOF(mu,a_f,e_f,th0_f, thf_f);

% total time lapse:
delta_t_tot = delta_t1+delta_t2+delta_t3+delta_t4;


%% ALTERNATIVE TRANSFER3 ~ REDUCED PLANE CHANGED Δv
clear;
clc;
close all;

% gravitational parameter:
mu = 398600; % [km^3/s^2]

% initial orbit data:
% Cartesian coordinates[km]
x = -5829.6395;
y = -6946.1419;
z = 4067.7268;
rvet = [x,y,z]';

% velocity vector components[km/s]
vx = 1.9470;
vy = -4.5900;
vz = -3.8490;
vvet = [vx,vy,vz]';

% Keplerian orbital parameters of initial orbit:
[a_i,e_i,i_i,OM_i,om_i,theta_i] = rv2parorb(rvet, vvet, mu);

% semi-latus rectum, periapsis and apoapsis of initial orbit:
p_i = a_i*(1-e_i^2); %[km]
r_p_i = a_i*(1-e_i); %[km]
r_a_i = a_i*(1+e_i); %[km]

% final orbit data:

% Keplerian orbital parameters of final orbit:
a_f = 13790.0000; %[km]
e_f = 0.2719; %[-]
i_f = 1.0470; % [rad]
OM_f = 2.5590; % [rad]
om_f = 0.4614; % [rad]
theta_f = 0.2578; % [rad]

% semi-latus rectum, periapsis and apoapsis of final orbit:
p_f = a_f*(1-e_f^2); %[km]
r_p_f = a_f*(1-e_f); %[km]
r_a_f = a_f*(1+e_f); %[km]

r_star_vet = []; %vector that will be filled with star radiuses 
delta_v_tot = []; %vector that will be filled with total velocity difference
% that belong to each star radius

% cycle to find optimal star radius
for r_star = r_a_i:1:100000
    r_star_vet = [r_star_vet, r_star];
    % Keplerian parameters for star orbit:
    i_star = i_i; %[rad]
    OM_star = OM_i; %[rad]
    om_star = om_i; %[rad]
    a_star = (r_star+r_p_i)/2; %[km]
    e_star = (r_star-r_p_i)/(r_star+r_p_i); %[-]
    
    % semi-latus rectum or star orbit:
    p_star = a_star*(1-e_star^2); %[km]
    
    % calculate velocity at pericenter respectively in initial orbit and
    % star orbit:
    v_transv_p_i = sqrt(mu/p_i)*(1+e_i); %[km/s]
    v_transv_p_star = sqrt(mu/p_star)*(1+e_star); %[km/s]
    % velocity difference for tangent manoeuvre to move from initial to
    % star orbit:
    deltav_1 = v_transv_p_star - v_transv_p_i; % [km/s]

    % true anomaly where we circularize star orbit:
    theta_circ = pi;
    [deltav_circ, r_circ,e_circ, i_circ, OM_circ, om_circ, theta_i_circ] = ellipse2circ(mu,a_star,e_star,i_star,OM_star,theta_circ);

    % plane change manoeuvre:
    [deltav_cp,om_cp,theta_cp1, theta_cp2] = changeOrbitalPlane(mu,r_circ,e_circ,i_circ,OM_circ,om_circ,i_f,OM_f,theta_i_circ);

    % distinguish two cases:
    if r_circ > r_a_f
        % "try orbit" that has r_circ as apoapsis and r_a_f as periapsis:
        a_try = (r_circ+r_a_f)/2;
        e_try = (r_circ-r_a_f)/(r_circ+r_a_f);
        p_try = a_try*(1-e_try^2);

        % velocity of circular orbit and velocity at try orbit apocenter:
        v_transv_circ = sqrt(mu/r_circ);
        v_transv_a_try = sqrt(mu/p_try)*(1-e_try);
        
        % velocity difference to move from circular orbit to try orbit:
        deltav_3 = v_transv_a_try -v_transv_circ;

        % velocity at try orbit pericenter and at final orbit apocenter:
        v_transv_p_try = sqrt(mu/p_try)*(1+e_try);
        v_transv_a_f = sqrt(mu/p_f)*(1-e_f);

        % velocity difference to move from try orbit into the final orbit:
        deltav_4 = v_transv_a_f - v_transv_p_try;

    elseif r_circ < r_a_f
        % "try orbit" that has r_a_f as apoapsis and r_circ as periapsis:
        a_try = (r_circ+r_a_f)/2;
        e_try = (r_a_f-r_circ)/(r_circ+r_a_f);
        p_try = a_try*(1-e_try^2);

        % velocity of circular orbit and velocity at try orbit pericenter:
        v_transv_circ = sqrt(mu/r_circ);
        v_transv_p_try = sqrt(mu/p_try)*(1+e_try);

        % velocity difference to move from circular orbit to try orbit:
        deltav_3 = v_transv_p_try -v_transv_circ;

        % velocity at try orbit apocenter and at final orbit apocenter:
        v_transv_a_try = sqrt(mu/p_try)*(1-e_try);
        v_transv_a_f = sqrt(mu/p_f)*(1-e_f);

        % velocity difference to move from try orbit into the final orbit:
        deltav_4 = v_transv_a_f - v_transv_a_try;
    end

    % total velocity difference:
    delta_v_tot = [delta_v_tot, abs(deltav_1)+abs(deltav_circ)+abs(deltav_cp)+abs(deltav_3)+abs(deltav_4)];
end

plot(r_star_vet,delta_v_tot);
hold on
grid on
xlabel('star radius [km]');
ylabel('total velocity difference [km/s]');

% after more refined discretisation we can assure that optimal radius is
% equal to final orbit apoapsis:
r_opt = r_a_f;

i_star = i_i;
OM_star = OM_i;
om_star = om_i;
a_star = (r_opt+r_p_i)/2;
e_star = (r_opt-r_p_i)/(r_opt+r_p_i);

p_star = a_star*(1-e_star^2);

v_transv_p_i = sqrt(mu/p_i)*(1+e_i);
v_transv_p_star = sqrt(mu/p_star)*(1+e_star);

deltav_1 = v_transv_p_star - v_transv_p_i;

theta_circ = pi;
[deltav_circ, r_circ, e_circ, i_circ, OM_circ, om_circ, theta_i_circ] = ellipse2circ(mu,a_star,e_star,i_star,OM_star,theta_circ);

[deltav_cp,om_cp,theta_cp1, theta_cp2] = changeOrbitalPlane(mu,r_circ,e_circ,i_circ,OM_circ,om_circ,i_f,OM_f,theta_i_circ);

v_circ = sqrt(mu/r_circ);
v_transv_a_f = sqrt(mu/p_f)*(1-e_f);

deltav_4 = v_transv_a_f - v_circ;

% since r_opt is equal to r_a_f we don't have a "try orbit" so we have not
% a deltav_3 as in the cycle
delta_v_tot = abs(deltav_1)+abs(deltav_circ)+abs(deltav_cp)+abs(deltav_4);

% Graphical representation:
figure
% declare initial,final true anomaly and discretisation step for initial
% orbit:
th0_i = theta_i;
thf_i = 2*pi;
dth_i = 0.01;
% initial orbit plot:
color = 1;
plotOrbit(a_i, e_i, i_i, OM_i, om_i, th0_i, thf_i, dth_i, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for star
% orbit:
th0_star = 0;
thf_star = pi;
dth_star = 0.01;
% star orbit plot:
color = 2;
plotOrbit(a_star, e_star, i_star, OM_star, om_star, th0_star, thf_star, dth_star, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for circular
% orbit:
th0_circ = pi+om_i;
thf_circ = theta_cp1;
dth_circ = 0.01;
% circular orbit plot:
color = 3;
plotOrbit(r_circ, e_circ, i_circ, OM_circ, om_circ, th0_circ, thf_circ, dth_circ, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for plane
% changed circular orbit:
th0_cp = theta_cp1 + om_cp - 2*pi;
thf_cp = 3*pi + (om_f);
dth_cp = 0.01;
% plane changed circular orbit plot:
color = 4;
plotOrbit(r_circ, e_circ, i_f, OM_f, om_circ, th0_cp, thf_cp, dth_cp, mu,color);
hold on

% declare initial,final true anomaly and discretisation step for final
% orbit:
th0_f = pi;
thf_f = 2*pi+theta_f;
dth_f = 0.01;
% final orbit plot:
color = 6;
plotOrbit(a_f, e_f, i_f, OM_f, om_f, th0_f, thf_f, dth_f, mu,color);
hold on

% earth plot:
earth3D

legend('initial orbit','initial orbit pericenter','star orbit','circularization point','circular orbit','plane change point','circular plane changed orbit','tangent manoeuvre point','final orbit','target point')

% calculate time lapses:

% time lapse from relapse point to tangent manoeuvre point:
[delta_t1] = TOF(mu, a_i, e_i, th0_i, thf_i);

% time lapse from tangent manoeuvre point to circularization point:
[delta_t2] = TOF(mu, a_star, e_star, th0_star, thf_star);

% time lapse from circularization point to plane change point:
[delta_t3] = TOF(mu, r_circ, e_circ, th0_circ, thf_circ);

% time lapse from plane change point to tangent manoeuvre point to change
% shape into final elliptical orbit:
[delta_t4] = TOF(mu, r_circ, e_circ, th0_cp, thf_cp);

% time lapse from tangent manoeuvre to final target point:
[delta_t5] = TOF(mu, a_f, e_f, th0_f, thf_f);

% total time lapse:
delta_t_tot = delta_t1+delta_t2+delta_t3+delta_t4+delta_t5;