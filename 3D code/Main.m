%% Lattice Boltzmann LBGK Test Code
% Adrian Harwood, The University of Manchester, UK
% Last Updated: March 2016

%% Initialisation

% Assumes D3Q19 and 19th trajectory is rest particle
clear
set(0,'DefaultFigureWindowStyle','docked')
close all

% Visuals
out_every = 10;

% Time
T = 3000;        % Number of time steps

% Lattice dimensions
N = 15;     % Number of x lattice sites
M = 15;     % Number of y lattice sites
K = 15;     % Number of z lattice sites

% Reynolds number
Re = 5;
lref = M;   % Reference length (in LUs) for Reynold number

% Lattice site coordinates dx = 1
xl = 1:1:N; yl = 1:1:M; zl = 1:1:K;
[x,y,z] = meshgrid(xl,yl,zl);   % Grid of coordinates for plotting

% Discrete velocities for D2Q9
c = [   1,	-1,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  0;...
        0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1,  0,  0,  0,  0,  0;...
        0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1,  0   ];
    
% Weights for D2Q9
w = [   1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,...
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,...
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,...
		1.0/3.0     ];

% Initialise matrices ((19 x) K x M x N)
feq = zeros(19,K,M,N);      % Equilibrium function
u = zeros(3,K,M,N);         % Velocity field
rho = zeros(K,M,N);         % Density field

% Forcing
forces_xy = zeros(3,K,M,N);     % Force field (Cartesian)
forces_i = zeros(19,K,M,N);     % Force field (Lattice)
F_extern = 0.0001;           % External force (could be gravity)
F_dir = 'x';                % String indicating the direction of the 
                            % external force ('x' or 'y' or 'z')

% Create typing matrix 'f' for fluid site 'b' for boundary
LatTyp = cell(K,M,N);
LatTyp(:,:,:) = {'f'};

% Initialise macroscopic quanitities
% Velocity Field
u0 = [0 .03 0];       % Initial Velocity Magnitude (x,y,z) (t = 0)
uref = .03;         % Reference velocity for Reynolds number
u(1,:,:,:) = u0(1);   % x-velocity field
u(2,:,:,:) = u0(2);   % y-velocity field
u(3,:,:,:) = u0(3);   % z-velocity field

% Density
rho0 = 1;               % Initial Density (arbitrary)
rho(:,:,:) = rho0;      % Uniform Field (arbitrary)

% Other quantities
cs = 1/sqrt(3);             % Sound speed on D3Q19 lattice
nu = uref * lref / Re;      % Kinematic viscosity from Reynolds
omega = 1 / ( (nu/cs^2) + .5 ); % Relaxation frequency related to viscosity
                                % for unit time step

% Initialise mesoscopic quantities
feq = Equilibrium(rho,w,c,u,cs,feq,N,M,K);
f = feq;
                                
% Information for user
disp(['Relaxation Time = ' num2str(1/omega)])
disp(['Kinematic Viscosity = ' num2str(nu)])
disp(['Reynolds Number = ' num2str(Re)])
disp('Running...')

% Visualise initial state
Visuals(N,M,K,x,y,z,u,rho,0);


%% LBM loop
tic
for t = 1 : T
    
    % Generate force vectors
    [forces_i, forces_xy] = Force(forces_i,forces_xy,LatTyp,...
        N,M,K,F_extern,F_dir,rho,u,omega,w,cs,c);
    
    % Compute collision
    [f, feq] = Collide(rho,w,c,u,cs,feq,N,M,K,f,omega,forces_i);

    % Apply boundary conditions
    f = Boundary(f,N,M,K,c,LatTyp);
    
    % Stream populations
    f = Stream(f,N,M,K,c);

    % Find macroscopic quantities
    [rho,u,ke] = Macroscopic(f,c,forces_xy,LatTyp);
    
    % Call visualisation
    if (mod(t,out_every) == 0)
        Visuals(N,M,K,x,y,z,u,rho,t);
    end
    
end
toc