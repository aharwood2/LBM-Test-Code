function [rho, u, ke] = Macroscopic(f,c)
% Compute macroscopic quantities

rho = sum(f,3);

for k = 1:9
    fu(:,:,k) = c(1,k) * f(:,:,k);
    fv(:,:,k) = c(2,k) * f(:,:,k);
end

u(:,:,1) = sum(fu,3) ./ rho;
u(:,:,2) = sum(fv,3) ./ rho;

% Total kinetic energy
ke = 0.5*sum(sum(sum(sqrt(u(:,:,1).^2 + u(:,:,2).^2) )));


end