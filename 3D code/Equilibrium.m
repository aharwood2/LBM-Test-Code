function [feq] = Equilibrium(rho,w,c,u,cs,feq,N,M,K)
% Loop through the lattice points to compute equilibrium based on
%       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4

for i = 1:N
    for j = 1:M
        for k = 1:K
        
            for v = 1:19
                
                % Compute c_ia * u_a which is actually the dot product of c
                % and u
                A = dot( c(:,v), u(:,k,j,i) );

                % Compute Q_iab u_a u_b = 
                %          (c_a^2 - cs^2)u_a^2 + 2c_a c_b u_a u_b + 
                %                   (c_b^2 - cs^2)u_b^2
                B = (c(1,v)^2 - cs^2) * u(1,k,j,i)^2 + ...
                    (c(2,v)^2 - cs^2) * u(2,k,j,i)^2 + ...
                    (c(3,v)^2 - cs^2) * u(3,k,j,i)^2 + ...
                    2*c(1,v)*c(2,v)*u(1,k,j,i)*u(2,k,j,i) + ...
                    2*c(1,v)*c(3,v)*u(1,k,j,i)*u(3,k,j,i) + ...
                    2*c(2,v)*c(3,v)*u(2,k,j,i)*u(3,k,j,i);

                % Compute f^eq
                feq(v,k,j,i) = rho(k,j,i) * w(v) * ...
                    (1 + (A / cs^2) + (B / (2*cs^4)));

            end
        
        end
    end
end
    
end
    
