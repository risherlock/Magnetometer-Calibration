function [] = plot_sphere(centre,radius)
% Plots sphere 
%
% Input:
%   centre = A vector [xc,yc,zc] representating centre of sphere. 
%   radius = Radius of sphere.
%
% 2020/06/04

x_c = centre(1);
y_c = centre(2);
z_c = centre(3);

% Circle parametric equation
theta       = linspace(0,2*pi,20);
phi         = linspace(0,pi,20);
[theta,phi] = meshgrid(theta,phi);
rho         = radius;

x = x_c + rho*sin(phi).*cos(theta);
y = y_c + rho*sin(phi).*sin(theta);
z = z_c + rho*cos(phi);

mesh(x,y,z,'FaceAlpha','0.4')
axis equal
hold on;
end
