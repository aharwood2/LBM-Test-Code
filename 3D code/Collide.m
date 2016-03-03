function [fnew, feq] = Collide(rho,w,c,u,cs,feq,N,M,K,f,omega,force_i)
% Compute the new ditribution functions using LBKG SRT collision

% Update equilibrium
feq = Equilibrium(rho,w,c,u,cs,feq,N,M,K);

% Recompute distribution function f
fnew = (-omega * (f - feq) ) + f + force_i;
    
end
    
