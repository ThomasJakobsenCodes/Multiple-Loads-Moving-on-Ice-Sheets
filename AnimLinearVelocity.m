% Computing the deflection of the surface of an ice sheet due to a load
% moving with constant linear velocity for consecutive time values
clear
clc
run('Constants.m'); % Loads the constant values

% Variables
% Spatial grid
N = 500; % Number of grid points
L = 2000; % Length of grid
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
time_end = 20;
time_steps = 100;
time = linspace(0, time_end, time_steps);
T = length(time);

% Velocity
velocity = [40, 10];

% Load mass
mass = 2000; % mass of the load in kg

% Computing the displacement
eta_hat = zeros(N, N, length(time));
for t = 1:T
    eta_hat(:, :, t) = compute_eta_hat(Xi_x, Xi_y, X, Y, Xi_norm, time(t), velocity, mass);
end
eta = zeros(N, N, length(time));
for t = 1:T
    eta(:, :, t) = real(ifft2((eta_hat(:, :, t))));
end

% Plotting
for t = 1:T
    clf
    surf(X, Y, eta(:, :, t))
    axis([min(X,[],"all") max(X,[],"all") min(Y,[],"all") max(Y,[],"all") min(eta,[],"all") max(eta,[],"all")]);
    title("Surface Displacement at Time ="+time(t))
    xlabel('X-axis (m)');
    ylabel('Y-axis (m)');
    zlabel('Displacement (m)');
    colorbar;
    shading interp;
    drawnow
end

% Function definitions
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

function result = integral(a, velocity_dot_xi, t)
    result = (exp(-1i * t * velocity_dot_xi) - exp(-t * a)) ...
             ./ (a - 1i * velocity_dot_xi);
end


function result = compute_eta_hat(xi_x, xi_y, x, y, xi_norm, t, velocity, mass)
    g_0 = compute_g_0(xi_norm);
    k = compute_k(xi_norm);
    r = compute_r(g_0, k);
    u = compute_u(xi_norm, g_0, k, r);

    intial_w_hat = compute_w_hat(mass, k, x, y);
    velocity_dot_xi = velocity(1) * xi_x + ...
                      velocity(2) * xi_y;
    integral1 = integral(r - 1i * u, velocity_dot_xi, t);
    integral2 = integral(r + 1i * u, velocity_dot_xi, t);
    result = 1i * g_0 ./ (2 * u) .* intial_w_hat .* (integral1 - integral2);
end
