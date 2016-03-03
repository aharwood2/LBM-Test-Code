function [forces_i, forces_xy] = Force(forces_i, forces_xy, LatTyp, ...
    N, M, ext_f, f_dir, rho, u, omega, w, cs, c)
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
forces_i(:,:,:) = 0;

% Reset Cartesian force vector on every grid site
forces_xy(:,:,:) = 0;

% Else, compute forces

% Get direction from string
if (strcmp(f_dir,'x'))
    dir = 1;
elseif (strcmp(f_dir,'y'))
    dir = 2;
else
    % Not a valid direction
    ext_f = 0;
    dir = 1;
end

% Loop over grid and add forces for each direction
for i = 1:N
    for j = 1:M
        
        % Only apply to non-solid sites
        if (~strcmp(LatTyp(j,i),'b'))

            % Add gravity to force vector
            forces_xy(j,i,dir) = forces_xy(j,i,dir) + ...
                (rho(j,i) * ext_f);

            % Now compute force_i components from Cartesian force vector
            for v = 1 : 9                

                % Reset beta_v
                beta_v = 0.0;

                % Compute the lattice forces based on Guo's scheme
                lambda_v = (1 - 0.5 * omega) * ( w(v) / (cs^2) );

                % Dot product
                beta_v = beta_v + (c(1,v) * u(j,i,1) + ...
                    c(2,v) * u(j,i,2));
                beta_v = beta_v * (1 / (cs^2));

                % Compute force on lattice X and Y contributions
                forces_i(j,i,v) = forces_i(j,i,v) + ...
                        forces_xy(j,i,1) * ...
                        (c(1,v) * (1 + beta_v) - u(j,i,1));
                forces_i(j,i,v) = forces_i(j,i,v) + ...
                        forces_xy(j,i,2) * ...
                        (c(2,v) * (1 + beta_v) - u(j,i,2));

                % Multiply by lambda_v
                forces_i(j,i,v) = forces_i(j,i,v) * lambda_v;

            end

        end

    end

end

end