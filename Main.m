% Lattice Boltzmann LBGK Test Code
% Adrian Harwood
% Last Updated: January 2015

% Assumes D2Q9 and 9th trajectory is rest particle
clear all
set(0,'DefaultFigureWindowStyle','docked')
close all
tic

% Time
T = 2; % Number of time steps
dt = 1; % Time step size
tcount = 0; % Reset counter

% Lattice dimensions
N = 32;  % Number of x lattice sites
M = 16;  % Number of y lattice sites

% Lattice site coordinates dx = 1
xl = 1:1:N; yl = 1:1:M;
[x,y] = meshgrid(xl,yl); % Grid of coordinates for plotting

% Discrete velocities for D2Q9
c = [1 0;...
    1 1;...
    0 1;...
    -1 1;...
    -1 0;...
    -1 -1;...
    0 -1;...
    1 -1;...
    0 0]';

% Weights for D2Q9
w = [1/9 1/36 1/9 1/36 1/9 1/36 1/9 1/36 4/9];

% Initialise population matrices (N x M (x 9))
f = zeros(M,N,9);       % Distribution function
feq = zeros(M,N,9);     % Equilibrium function
u = zeros(M,N,2);       % Velocity field
rho = zeros(M,N);       % Density field

% Create typing matrix 'f' for Fluid site
LatTyp = cell(M,N);
LatTyp(:,:) = {'f'};
% Boundary labels
LatTyp(1,:) = {'b'}; LatTyp(M,:) = {'b'};
LatTyp(2:M-1,1) = {'b'}; LatTyp(2:M-1,N) = {'b'};

% Initialise macroscopic quanitities
% Velocity Field
u0 = [.02 0];       % Initial Velocity Magnitude (x,y) (t = 0)
u(:,:,1) = u0(1); % x-velocity field
u(:,:,2) = u0(2); % y-velocity field

% Density
rho0 = 1;           % Initial Density (arbitrary)
rho(:,:) = rho0;    % Uniform Field (arbitrary)

% Viscosity
nu = .02;           % Kinematic viscosity

% Other quantities
cs = 1/sqrt(3);   % Sound speed on D2Q9 lattice
omega = 1 / ( (nu/cs^2) + .5 ); % Relaxation frequency related to viscosity
                                % for unit time step
                                
disp(['Relaxation Time = ' num2str(1/omega)])
disp(['Kinematic Viscosity = ' num2str(nu)])

% Reynolds number (domain width as length scale)
L = N;
Re = norm(u0) * L / nu;

% LBM
for t = dt:dt:T
    
    % Compute collision
    [f, feq] = Collide(rho,w,c,u,cs,feq,N,M,f,omega);

    % Apply boundary conditions
    f = Boundary(f,u,u0,rho,N,M);
    
    % Stream populations
    f = Stream(f,N,M,c);

    % Find macroscopic quantities
    [rho,u,ke] = Macroscopic(f,c);
    
    % Increment counter
    tcount = tcount + 1;
    
    % Call visualisation
    Visuals(x,y,u)
    
end
toc