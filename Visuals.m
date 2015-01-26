function Visuals(x,y,u)
% Function to display quiver plot of velocity

figure
quiver(x,y,u(:,:,1),u(:,:,2))
axis equal
axis tight