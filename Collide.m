function [fnew, feq] = Collide(rho,w,c,u,cs,feq,N,M,f,omega)
% Loop through the lattice points to compute the new ditribution functions.
% Equilibrium based on:
%       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4

for i = 1:N
    for j = 1:M
        
        for k = 1:9
            % Compute c_ia * u_a which is actually the dot product of c and
            % u
            A = dot(c(:,k),reshape(u(j,i,:),2,1));

            % Compute Q_iab u_a u_b = 
            %          (c_a^2 - cs^2)u_a^2 + 2c_a c_b u_a u_b + 
            %                   (c_b^2 - cs^2)u_b^2
            B = (c(1,k)^2 - cs^2) * u(j,i,1)^2 + ...
                2*c(1,k)*c(2,k)*u(j,i,1)*u(j,i,2) + ...
                (c(2,k)^2 - cs^2) * u(j,i,2)^2;

            % Compute f^eq
            feq(j,i,k) = rho(j,i) * w(k) * ...
                (1 + (A / cs^2) + (B / (2*cs^4)));
            
            % Alternative form for feq I have seen published
%             feq2(j,i,k) = rho(j,i) * w(k) * ...
%                 ( 1 + ...
%                 ((3/2)*dot(c(:,k),reshape(u(j,i,:),2,1))) + ...
%                 ((9/2)*(dot(c(:,k),reshape(u(j,i,:),2,1))^2)) - ...
%                 ((3/2)*norm(reshape(u(j,i,:),2,1))^2) );
            
        end
        
    end
end

% Recompute distribution function f
fnew = (-omega * (f - feq) ) + f;
    
end
    
