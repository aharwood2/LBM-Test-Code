function [rho, u, ke] = Macroscopic(f,c,force_xy)
% Compute macroscopic quantities
fu = zeros(size(f,1),size(f,2));
fv = zeros(size(f,1),size(f,2));

% Density
rho = sum(f,3);

% Mass flux
for k = 1:9
    fu = fu + (c(1,k) * f(:,:,k));
    fv = fv + (c(2,k) * f(:,:,k));
end

% Add forces
fu = fu + (0.5 * rho .* force_xy(:,:,1));
fv = fv + (0.5 * rho .* force_xy(:,:,2));

% Get velocity
u(:,:,1) = fu ./ rho;
u(:,:,2) = fv ./ rho;

% Total kinetic energy
ke = 0.5*sum(sum(sum(sqrt(u(:,:,1).^2 + u(:,:,2).^2) )));


end