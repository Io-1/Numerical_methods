u0=@(x) x.^2+2;
g=@(t) 3*t+2;
l=@(t) 3;
sourc=@(x,t) x+t;

a=2;
mu=1/6;
h=0.025;
r=mu*h^2/a^2

U=Theta_method([0,1],h,[0,1.5],r,a,u0,g,l,sourc,0);




