function Visuals(N,M,K,x,y,z,u,rho,t)
% Function to display stuff

figure

% Just Quiver
quiver3(x,y,z,...
    reshape(u(1,:,:,:),size(x)),...
    reshape(u(2,:,:,:),size(x)),...
    reshape(u(3,:,:,:),size(x)));
axis equal
axis tight
title(['Time Step  = ' num2str(t)])



% subplot(3,1,1), quiver3(x,y,z,...
%     reshape(u(1,:,:,:),size(x)),...
%     reshape(u(2,:,:,:),size(x)),...
%     reshape(u(3,:,:,:),size(x)));
% axis equal
% axis tight
% title(['Time Step  = ' num2str(t)])
% 
% umag = sqrt(u(1,:,:,:).^2 + u(2,:,:,:).^2 + u(3,:,:,:).^2);
% subplot(3,1,2), scatter3(reshape(x,N*M*K,1),reshape(y,N*M*K,1),reshape(z,N*M*K,1),...
%                 1,reshape(umag,N*M*K,1));
% colorbar
% axis([-0.5 N + 0.5 -0.5 M + 0.5 -0.5 K + 0.5])
% axis equal
% title(['Time Step  = ' num2str(t)])
% 
% subplot(3,1,3), scatter3(reshape(x,N*M*K,1),reshape(y,N*M*K,1),reshape(z,N*M*K,1),...
%                 1,reshape(rho,N*M*K,1));
% colorbar
% axis([-0.5 N + 0.5 -0.5 M + 0.5 -0.5 K + 0.5])
% axis equal
% title(['Time Step  = ' num2str(t)])
drawnow

end