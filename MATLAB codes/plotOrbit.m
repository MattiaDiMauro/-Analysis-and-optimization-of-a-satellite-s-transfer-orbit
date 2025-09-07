function plotOrbit(a, e, i, OM, om, th0, thf, dth, mu,color)

% Function to plot a 3D orbit in the intertial frame:
%
%
% plotOrbit(a, e, i, OM, om, th0, thf, dth, mu, color)
%
% Input arguments:
% -------------------------------------------------------------------------
% a         [1x1]   semi-major axis                             [km]
% e         [1x1]   eccentricity                                [-]
% i         [1x1]   inclination                                 [rad]
% OM        [1x1]   RAAN                                        [rad]
% om        [1x1]   pericenter anomaly                          [rad]
% th0       [1x1]   initial true anomaly                        [rad]
% thf       [1x1]   final true anomaly                          [rad]
% dth       [1x1]   true anomaly discretisation step            [rad]
% mu        [1x1]   gravitational parameter                     [km^3/s^2]
% color     [1x1]   color of the orbit line                     [-]
% -------------------------------------------------------------------------

if color == 1
    c = [0 0.4470 0.7410];
elseif color == 2
    c = [0.8500 0.3250 0.0980];
elseif color == 3
    c = [0.9290 0.6940 0.1250];
elseif color == 4
    c = [0.4940 0.1840 0.5560];
elseif color == 5
    c = [0.4660 0.6740 0.1880];
elseif color == 6
    c = [0.3010 0.7450 0.9330];
elseif color == 7
    c = [0.6350 0.0780 0.1840];
end

rplot = [];

for j = th0:dth:thf
    theta = j;
    [rr,vv] = parorb2rv( a,e,i,OM,om,theta,mu);
    rplot = [rplot, rr];
end

plot3(rplot(1,:),rplot(2,:),rplot(3,:), 'Color', c);
hold on
grid on
hold on
% add a pointmarker in orbital point corresponding to the thf true
% anomaly
plot3(rplot(1,end),rplot(2,end),rplot(3,end),'.','MarkerSize',13,'MarkerFaceColor',c,'MarkerEdgeColor',c);
xlabel('[km]');
ylabel('[km]');
zlabel('[km]');
end