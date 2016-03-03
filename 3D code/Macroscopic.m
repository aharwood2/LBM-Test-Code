function [rho, u, ke] = Macroscopic(f,c,force_xyz,LatTyp)
% Compute macroscopic quantities
fu = zeros(1,size(f,2),size(f,3),size(f,4));
fv = zeros(1,size(f,2),size(f,3),size(f,4));
fw = zeros(1,size(f,2),size(f,3),size(f,4));

% Density
rho(:,:,:) = sum(f,1);

% Mass flux
for v = 1:19
    fu = fu + (c(1,v) * f(v,:,:,:));
    fv = fv + (c(2,v) * f(v,:,:,:));
    fw = fw + (c(3,v) * f(v,:,:,:));
end

% Add forces
fu = fu + (0.5 * reshape(rho,[1 size(rho)]) .* force_xyz(1,:,:,:));
fv = fv + (0.5 * reshape(rho,[1 size(rho)]) .* force_xyz(2,:,:,:));
fw = fw + (0.5 * reshape(rho,[1 size(rho)]) .* force_xyz(3,:,:,:));

% Get velocity
u(1,:,:,:) = fu ./ reshape(rho,[1 size(rho)]);
u(2,:,:,:) = fv ./ reshape(rho,[1 size(rho)]);
u(3,:,:,:) = fw ./ reshape(rho,[1 size(rho)]);

% Apply resets
u( repmat( strcmp(LatTyp,{'b'}), 3,1,1,1 ) ) = 0;

% Total kinetic energy
ke = 0.5*sum(sum(sum(sum( sqrt(u(1,:,:,:).^2 + u(2,:,:,:).^2 + u(3,:,:,:).^2) ))));


end