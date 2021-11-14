%Modified and Improved by Mainak Mallick
clear all
close all
clc
cd0 = 0.0072;
%condition%
rho = 0.8194;
Vinf = 25;
alpha = 6*(pi/180);
%wing%
lambda = .4; %taper ratio
b_half =6.25;
Cr = .8; %root chord
a0 = 6.393657747; %lift slope
part = 100;
z = 0;
for i = 1 : (part/4)+2
    theta(i) = pi/2 - z*(pi/(2*part));
    Y(i) = b_half*cos(theta(i));
    nY(i) = Y(i)/b_half;
    c(i) = .8; %rectangular wing 
    %c(i) = Cr-((i-1)*Cr*(1-lambda)/(part));
    %c(i) = Cr*((1 - lambda) * (b_half - Y(i))/ (b_half + lambda)); %taper wing
    %c(i) = Cr*(1 - ((Y(i)/b_half)^2))^0.5; %elliptic wing
    vr(i) = Vinf;
    mu(i) = (a0*c(i)/(8*b_half));
    thetatwist(i) = 0; %no twist
    %thetatwist(i) = -2/b_half*Y(i)*pi/180; %geometric twist
    alphacl0(i) = -4.8*pi/180;
    alphaY(i) = thetatwist(i) + alpha;
    z = z + 1;
end
for i = (part/4)+3: part
    theta(i) = pi/2 - z*(pi/(2*(part)));
    
    Y(i) = b_half*cos(theta(i));
    nY(i) = Y(i)/b_half;
    %c(i) = 1; %rectangular wing 
    %c(i) = Cr-((i-1)*Cr*(1-lambda)/(part-((2*part/5)+1)));
    c(i) = Cr*((1 - lambda-.04) +((b_half-.85 - Y(i))/ (b_half + lambda))); %taper wing
    %c(i) = Cr*(1 - ((Y(i)/b_half)^2))^0.5; %elliptic wing
    vr(i) = Vinf;
    mu(i) = (a0*c(i)/(8*b_half));
    thetatwist(i) = 0; %no twist
    %thetatwist(i) = -2/b_half*Y(i)*pi/180; %geometric twist
    alphacl0(i) = -4.8*pi/180;
    alphaY(i) = thetatwist(i) + alpha;
    z = z + 1;
end
for i = 1 : part-1
    dY(i) = Y(i+1) - Y(i);
end
dY(part) = 0; 
n = 0;
for i = 1 : part
    for j = 1 : part
        a(i,j) = sin((2*n+1)*theta(i))*((2*n+1)*mu(i) + sin(theta(i)));
        n = n + 1;
    end
    n = 0;
end
area = 0;
for i = 1 : part 
    area = area + c(i)*Y(i);
end
for i = 1 : part
    b(i) = mu(i)*(alphaY(i) - alphacl0(i)) * sin(theta(i));
end
x = linsolve(a,b');
n = 0; xsum = 0;
for i = 1 : part
    for j = 1 : part
        xsum = xsum + 4*b_half*vr(i)*x(j)*sin((2*n+1)*theta(i));
        n = n + 1;
    end
    gamma(i) = xsum;
    lift(i) = gamma(i)*rho*vr(i); %kutta 
    cl(i) = lift(i)/(0.5*rho*Vinf^2*c(i)); %lift coefficient
    xsum = 0; n = 0;
end
xsum = 0; n = 0;
for i = 1 : part
    for j = 1: part
        xsum = xsum + (2*n+1) * x(j) * sin(theta(i) * (2*n+1))/sin(theta(i));
        n = n + 1;
    end
    alphai(i) = xsum;
    xsum = 0; n = 0;
end

for i = 1 : (2*part/5)
    drag(i) = -lift(i)*sin(alphai(i)) - 0.5*rho*vr(i)^2 * c(i)*dY(i)*cd0;
    wd(i) = (vr(i))*sin(alphai(i));
end
for i = (2*part/5)+1 : part
    drag(i) = -lift(i)*sin(alphai(i)) - 0.5*rho*vr(i)^2 * c(i)*dY(i)*cd0;
    wd(i) = vr(i)*sin(alphai(i));
end
cfv = 0; cft = 0; torq = 0;
for i = 1 : part
    cfv = cfv + lift(i)*dY(i);
    cft = cft * drag(i)*dY(i);
end
liftgrams = cfv*1000/9.8;
idraggrams = cft*1000/9.8;
LbyD = abs(cfv/cft);

figure(1)
plot(Y,gamma)
xlabel('span'); ylabel('gamma');

figure(2)
plot(Y,lift)
xlabel('span'); ylabel('cl/b');

figure(3)
plot(Y,-wd)
xlabel('span'); ylabel('downwash m/s');
figure(4)
plot(Y,c)
xlabel('span'); ylabel('chord');