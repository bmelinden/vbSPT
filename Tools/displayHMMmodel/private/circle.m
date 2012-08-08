
function circle(center,radius,nop)
THETA=linspace(0,2*pi,nop);
RHO=ones(1,nop)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
zData = ones(size(X));
patch(X,Y, zData);
axis equal;
end