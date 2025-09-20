function [ ex, ey, ez ] = ellipsoidrand (a,b,c,n)
n = floor(n);
%Normalize the ellipsoid to have a smallest axis of 1.
ai = [a,b,c]./min([a,b,c]);

%Generate points with uniform density on a sphere surface
u = rand(n,1); theta = 2*pi*u;
v = rand(n,1); phi = asin(2*v-1);
x = cos(theta).*cos(phi); y=sin(theta).*cos(phi); z=sin(phi);

%Select points with probability inversely related to how
%far they are from the surface of the original sphere.
%The acceptance rate appears to be about 50% in the worst
%case of an elongated cigar shape.
t = rand(size(u));
idx = find((x/ai(1)).^2 + (y/ai(2)).^2 + (z/ai(3)).^2 >= t.^2);
acceptance_rate = 100*size(idx)/size(t);

%stretch the selected points onto the ellipse
ex = x(idx)*a; ey=y(idx)*b; ez = z(idx)*c;
