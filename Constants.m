classdef Constants
    % These values are gathered from
    % T. Takizawa, “Deflection of a lfoating sea ice sheet induced by a moving load,” Cold Regions Science and Technologi, vol. 11, no. 2, pp. 171-180, 1985. 
    % K. Johnsen, “Wave response to an arbitrary motion of a load on an ice plate,” University of Bergen, 2022.
    properties (Constant = true)
        g = 9.81;           % gravity
        rho_water = 1026;   % water density
        rho_ice = 917;      % ice density
        h_water = 6.8;      % water depth
        h_ice = 0.17;       % ice thickness
        E = 5.1 * 10^8;     % Young's modulus
        poisson = 1 / 3;     % Poisson's ratio
        D = Constants.E * Constants.h_ice^3 / (12 * (1 - Constants.poisson^2)); % flexural rigidity
        kappa = Constants.D / (Constants.rho_water * Constants.g); % hydroelastic parameter
        b = 0.41 * 2 * sqrt(Constants.g * Constants.rho_water * ...
                            Constants.rho_ice * Constants.h_ice); % damping coefficient
    end
end