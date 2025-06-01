% Computing the deflection of the surface of an ice sheet due to multiple
% loads moving in an arbitrary path
clear
clc
run('Constants.m'); % Loads the constant values

% Variables
% Spatial grid
N = 500;   % Number of grid points
L = 6000; % Length of grid
x = linspace(-L / 2, L / 2, N);
y = linspace(-L / 2, L / 2, N);
[X, Y] = meshgrid(x, y);

% Fourier grid
xi_x = fftshift(2 * pi / L * (-floor(N / 2):ceil(N / 2) - 1));
xi_y = fftshift(2 * pi / L * (-floor(N / 2):ceil(N / 2) - 1));
[Xi_x, Xi_y] = meshgrid(xi_x, xi_y);
Xi_norm = sqrt(Xi_x.^2 + Xi_y.^2);
epsilon = 10^-12;
Xi_norm = max(Xi_norm, epsilon); % prevent division by 0 error

% Time
t = 75;

% Load mass
mass = 2000; % mass of the load in kg

% Computing the deflection
load_number = 10;
eta_s = zeros(N, N);
parfor load = 1:load_number
    load_progress = 1+load_number-load;
    disp("Computing load "+load_progress+" out of "+load_number)
    eta_hat = compute_eta_hat(Xi_x, Xi_y, X, Y, Xi_norm, t, load, mass);
    eta = real(ifft2(eta_hat));
    eta_s = eta_s + eta;
end

% Plotting
figure
surf(X, Y, eta_s)
title("Surface Displacement at Time ="+t)
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Displacement (m)');
colorbar;
shading interp;
end

% Function definitions
function result = path(t, load)
% PATH Parameterization of the path the loads moves in
%   result is a vector whose first component is the x-coordinate of the
%   load position and the second component is the y-coordinate. I.eg.
%   result=[x;y]
    v = 20;
    v_x_1 = 0;
    v_y_1 = v;
    v_x_2 = -v;
    v_y_2 = 0;
    s = 600; % length of the straight before the turn
    d = 100; % distance between carts
    t_standstill = 10;
    t_c = (s + (load - 1) * d) / v; % time the cart takes to get to the turn
    c_s = s + (load - 1) * d - v_y_1 * t_standstill; % distance from the cart to the turn
    rad = 300; % turn radius
    w = v / rad; % angular velocity
    T = pi / 2 / w; % period
    L_y = -c_s + v_y_1 * (t_c - t_standstill); % vertical endpoint of the first linear path
    C_y = L_y + rad; % vertical endpoint of the circular path

    if t < t_standstill
        result = [0; -c_s];
    elseif t < t_c
        result = [v_x_1 * (t - t_standstill); -c_s + v_y_1 * (t - t_standstill)];
    elseif t < (t_c + T)
        result = [-rad + rad * cos(w * (t - t_c)); L_y + rad * sin(w * (t - t_c))];
    else
        result = [-rad + v_x_2 * (t - (t_c + T)); C_y + v_y_2 * (t - (t_c + T))];
    end
end

function result = load_pressure(mass, x, y)
    sigma = 10;
    result = mass * 1 / (2 * pi * sigma) * exp(-(x.^2 + y.^2) / (2 * sigma^2));
end

function result = load_pressure_hat(mass, x, y)
    result = fft2(load_pressure(mass, x, y));  
end

function result = compute_w_hat(mass, k, x, y)
    result = load_pressure_hat(mass, x, y) ./ (Constants.rho_water * k);
end

function result = compute_g_0(xi_norm)
    result = xi_norm .* tanh(Constants.h_water * xi_norm);
end

function result = compute_k(xi_norm)
    result = 1 + (Constants.rho_ice * Constants.h_ice ...
            / Constants.rho_water) .* compute_g_0(xi_norm) ...
             .* (1 + (Constants.h_ice^2 / 12) * xi_norm.^2);
end

function result = compute_r(g_0, k)
    result = Constants.b / 2 / Constants.rho_water .* g_0 ./ k;
end

function result = compute_u(xi_norm, g_0, k, r)
    result = sqrt(Constants.g * (1 + Constants.kappa * xi_norm.^4) ...
                  .* g_0 ./ k - r.^2);
end

function result = compute_eta_hat(xi_x, xi_y, x, y, xi_norm, t, load, mass)
    g_0 = compute_g_0(xi_norm);
    k = compute_k(xi_norm);
    r = compute_r(g_0, k);
    u = compute_u(xi_norm, g_0, k, r);

    intial_w_hat = compute_w_hat(mass, k, x, y);
    a = r - 1i*u;
    b = r + 1i*u;
    integrand = @(tau) exp(-1 * (t - tau) * a - 1i * ...
                       ([1, 0] * path(tau, load) * xi_x + [0, 1] * path(tau, load) * xi_y)) ...
                       - exp(-1 * (t - tau) * b - 1i * ...
                       ([1, 0] * path(tau, load) * xi_x + [0, 1] * path(tau, load) * xi_y));
    compute_integral = integral(integrand, 0, t, 'ArrayValued', true);
    result = 1i * g_0 ./ (2 * u) .* intial_w_hat .* compute_integral;
end