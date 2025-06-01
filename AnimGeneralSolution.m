% Computing the deflection of the surface of an ice sheet due to a load
% moving in an arbitrary path for consecutive time values
clear
clc
run('Constants.m'); % Loads the constant values

% Variables
% Spatial grid
N = 500; % Number of grid points
L = 3000; % Length of grid
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
time_end = 100;
time_steps = time_end*2;
time = linspace(0, time_end, time_steps);
T = length(time);

% Load mass
mass = 2000; % mass of the load in kg

% Computing the displacement
eta_hat = zeros(N, N, length(time));
for t = 1:T
    eta_hat(:, :, t) = compute_eta_hat(Xi_x, Xi_y, X, Y, Xi_norm, time(t), time_end, mass);
    disp("Progress at "+(t/T*100)+"%")
end
eta = zeros(N, N, length(time));
for t = 1:T
    eta(:, :, t) = real(ifft2((eta_hat(:, :, t))));
end

% Saving a video of the plot as this code will take a long time to run
v = VideoWriter('myplot.avi');
v.FrameRate = 1;
open(v)

% Plotting
for t = 1:T
    surf(X, Y, eta(:, :, t))
    axis([min(X,[],"all") max(X,[],"all") min(Y,[],"all") max(Y,[],"all") min(eta,[],"all") max(eta,[],"all")]);
    title("Surface Displacement at Time ="+time(t))
    xlabel('X-axis (m)');
    ylabel('Y-axis (m)');
    zlabel('Displacement (m)');
    colorbar;
    shading interp;
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v)

% Function definitions
function result = path(t, time_end)
% PATH Parameterization of the path the loads moves in
%   result is a vector whose first component is the x-coordinate of the
%   load position and the second component is the y-coordinate. I.eg.
%   result=[x;y]
    T_half  = (time_end / 2 - time_end / 8);
    v_x = 0;
    v_y = 10;
    x_Linear_endpoint = v_x * T_half;
    y_Linear_endpoint = v_y * T_half;
    A = 300;
    if (t < time_end / 8)
        result = [0; 0];
    elseif (t < time_end / 2)
        result = [v_x * (t - time_end / 8); v_y * (t - time_end / 8)];
    else
        result = [x_Linear_endpoint - A + A * cos(pi/time_end * (t - time_end / 2)); y_Linear_endpoint + A * sin(pi/time_end * (t - time_end / 2))];
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

function result = compute_eta_hat(xi_x, xi_y, x, y, xi_norm, t, time_end, mass)
    g_0 = compute_g_0(xi_norm);
    k = compute_k(xi_norm);
    r = compute_r(g_0, k);
    u = compute_u(xi_norm, g_0, k, r);

    intial_w_hat = compute_w_hat(mass, k, x, y);
    a = r - 1i*u;
    b = r + 1i*u;
    integrand = @(tau) exp(-1 * (t - tau) * a - 1i * ...
                       ([1, 0] * path(tau, time_end) * xi_x + [0, 1] * path(tau, time_end) * xi_y)) ...
                       - exp(-1 * (t - tau) * b - 1i * ...
                       ([1, 0] * path(tau, time_end) * xi_x + [0, 1] * path(tau, time_end) * xi_y));
    compute_integral = integral(integrand, 0, t, 'ArrayValued', true);
    result = 1i * g_0 ./ (2 * u) .* intial_w_hat .* compute_integral;
end
