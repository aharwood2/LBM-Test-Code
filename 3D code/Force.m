function [forces_i, forces_xyz] = Force(forces_i, forces_xyz, LatTyp, ...
    N, M, K, ext_f, f_dir, rho, u, omega, w, cs, c)
% This routine computes the forces applied along each direction on the 
% lattice from Guo's 2002 scheme. The basic LBM must be modified in two 
% ways:
%       1)  The forces are added to the populations produced by the LBGK 
%           collision in the collision routine;
%       2)  dt/2 * F is added to the momentum in the macroscopic 
%           calculation. This has been done already in the other routines.
%
% The forces along each direction are computed according to:
%
%   F_i = (1 - 1/2*tau) * (w_i / cs^2) * 
%               ( (c_i - v) + (1/cs^2) * (c_i . v) * c_i ) . F
%
% where
%   F_i = force applied to lattice in i-th direction
%   tau = relaxation time
%   w_i = weight in i-th direction
%   cs = lattice sound speed
%   c_i = vector of lattice speeds for i-th direction
%   v = macroscopic velocity vector
%   F = Cartesian force vector
%
% In the following, we assume
%	lambda_i = (1 - 1/2*tau) * (w_i / cs^2)
%	beta_i = (1/cs^2) * (c_i . v)
%

% Reset lattice force vectors on every grid site
forces_i(:,:,:,:) = 0;

% Reset Cartesian force vector on every grid site
forces_xyz(:,:,:,:) = 0;

% Get direction from string
if (strcmp(f_dir,'x'))
    dir = 1;
elseif (strcmp(f_dir,'y'))
    dir = 2;
elseif (strcmp(f_dir,'z'))
    dir = 3;
else
    % Not a valid direction
    ext_f = 0;
    dir = 1;
end

% Loop over grid and add forces for each direction
for i = 1:N
    for j = 1:M
        for k = 1:K
        
            % Only apply to non-solid sites
            if (~strcmp(LatTyp(k,j,i),'b'))

                % Add gravity to force vector
                forces_xyz(dir,k,j,i) = forces_xyz(dir,k,j,i) + ...
                    (rho(k,j,i) * ext_f);

                % Now compute force_i components from Cartesian force vector
                for v = 1 : 19                

                    % Reset beta_v
                    beta_v = 0.0;

                    % Compute the lattice forces based on Guo's scheme
                    lambda_v = (1 - 0.5 * omega) * ( w(v) / (cs^2) );

                    % Dot product
                    beta_v = beta_v + (c(1,v) * u(1,k,j,i) + ...
                        c(2,v) * u(2,k,j,i) + c(3,v) * u(3,k,j,i));
                    beta_v = beta_v * (1 / (cs^2));

                    % Compute force on lattice X,Y,Z contributions
                    for d = 1 : 3
                        forces_i(v,k,j,i) = forces_i(v,k,j,i) + ...
                                forces_xyz(d,k,j,i) * ...
                                (c(d,v) * (1 + beta_v) - u(d,k,j,i));
                    end

                    % Multiply by lambda_v
                    forces_i(v,k,j,i) = forces_i(v,k,j,i) * lambda_v;

                end

            end

        end
    end
end

end