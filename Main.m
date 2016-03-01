%% Lattice Boltzmann LBGK Test Code
% Adrian Harwood, The University of Manchester, UK
% Last Updated: March 2016

%% Initialisation

% Assumes D2Q9 and 9th trajectory is rest particle
clear
set(0,'DefaultFigureWindowStyle','docked')
close all

% Visuals
c_scale = [0 1];        % Velocity scale for plotting

% Time
T = 100;        % Number of time steps

% Lattice dimensions
N = 50;  % Number of x lattice sites
M = 10;  % Number of y lattice sites

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

% Initialise matrices (N x M (x 9))
f = zeros(M,N,9);       % Distribution function
feq = zeros(M,N,9);     % Equilibrium function
u = zeros(M,N,2);       % Velocity field
rho = zeros(M,N);       % Density field

% Forcing
forces_xy = zeros(M,N,2);   % Force field (Cartesian)
forces_i = zeros(M,N,9);    % Force field (Lattice)
F_extern = 0.001;           % External force (could be gravity)
F_dir = 'x';                % String indicating the direction of the 
                            % external force (either 'x' or 'y')

% Create typing matrix 'f' for fluid site
LatTyp = cell(M,N);
LatTyp(:,:) = {'f'};
% Boundary labels 'b'
% LatTyp(1,:) = {'b'}; LatTyp(M,:) = {'b'};
% LatTyp(2:M-1,1) = {'b'}; LatTyp(2:M-1,N) = {'b'};

% Initialise macroscopic quanitities
% Velocity Field
u0 = [0 0];         % Initial Velocity Magnitude (x,y) (t = 0)
uref = .02;         % Reference velocity for Reynolds number
u(:,:,1) = u0(1);   % x-velocity field
u(:,:,2) = u0(2);   % y-velocity field

% Density
rho0 = 1;           % Initial Density (arbitrary)
rho(:,:) = rho0;    % Uniform Field (arbitrary)

% Reynolds number
Re = 10;
lref = M;   % Reference length (in LUs) for Reynold number

% Other quantities
cs = 1/sqrt(3);   % Sound speed on D2Q9 lattice
nu = uref * lref / Re;      % Kinematic viscosity from Reynolds
omega = 1 / ( (nu/cs^2) + .5 ); % Relaxation frequency related to viscosity
                                % for unit time step
                                
% Information for user
disp(['Relaxation Time = ' num2str(1/omega)])
disp(['Kinematic Viscosity = ' num2str(nu)])
disp(['Reynolds Number = ' num2str(Re)])
disp('Running...')


%% LBM loop
tic
for t = 1 : T
    
    % Generate force vectors
    [forces_i, forces_xy] = Force(forces_i,forces_xy,LatTyp,...
        N,M,F_extern,F_dir,rho,u,omega,w,cs,c);
    
    % Compute collision
    [f, feq] = Collide(rho,w,c,u,cs,feq,N,M,f,omega,forces_i);

    % Apply boundary conditions
    f = Boundary(f,u,u0,rho,N,M);
    
    % Stream populations
    f = Stream(f,N,M,c);

    % Find macroscopic quantities
    [rho,u,ke] = Macroscopic(f,c,forces_xy);
    
    % Call visualisation
    Visuals(x,y,u,t,c_scale);
    
end
toc