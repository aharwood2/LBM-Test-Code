function Visuals(x,y,u,t,c_scale)
% Function to display quiver plot of velocity

if t == 1
    figure
end
subplot(2,1,1), quiver(x,y,u(:,:,1),u(:,:,2))
axis equal
axis tight
title(['Time Step  = ' num2str(t)])

umag = sqrt(u(:,:,1).^2 + u(:,:,2).^2);
subplot(2,1,2), surf(x,y,umag,'EdgeColor','None')
view(2)
colorbar
caxis(c_scale);
axis equal
axis tight
title(['Time Step  = ' num2str(t)])
drawnow

end