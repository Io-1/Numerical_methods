function U= Theta_method(xspan,h,tspan,r,a,u0,g,l,f,TT)
%xspan is the span of space axis [x_0;L] as an array
%h is the step of space discretization as constant
%tspan is the span of space axis [0;T] as an array
%k is the step of time discretization as constant
%a is the constant
%u0 is the initial condition u(x,0)=u0(x) as a function handle
%g is the boundary condition u(x_0,t)=f(t) as a function handle
%l is the boundary condition u(L,t)=g(t) as a function handle
%f is the heat source
%TT is 0<=Theta<=1
%
%
clc

mu=a^2*r/h^2;

if TT<0 || TT>1
    disp('Theta is not in the bounds [0,1].')
    return
end

if 2*(1-TT)*mu <= 1
    disp('The maximum principle holds.')
else
    disp('The maximum principle does NOT hold.')
end

if TT>=0.5
    disp('Unconditionally stable (TT>=0.5).')
elseif 2*mu*(1-2*TT)<=1
    disp('L^2 stable solution.')
else
    disp('Unstable numerical solution.')
end

[X,T] = ndgrid(xspan(1):h:xspan(2), tspan(1):r:tspan(2));
N=length((xspan(1):h:xspan(2)));
M=length((tspan(1):r:tspan(2)));
U = ones(N,M);

U(:,1) = u0(X(:,1));
U(1,:) = g(T(1,:));
U(N,:) = l(T(N,:));

a(1)=0;
c(N-2)=0;
d=ones(1,N-2);

for k=1:M-1
    d(1)=(1-2*mu*(1-TT))*U(2,k)+mu*(1-TT)*(U(1,k)+U(3,k))+mu*TT*U(1,k+1)+(1-TT)*r*f(X(2),T(1,k))+TT*r*f(X(2),T(1,k+1));
    d(N-2)=((1-2*mu*(1-TT))*U(N-1,k)+mu*(1-TT)*(U(N-2,k)+U(N,k)))+mu*TT*U(N,k+1)+(1-TT)*r*f(X(N-1),T(1,k))+TT*r*f(X(N-1),T(1,k+1));
    a(2:N-2)=-mu*TT;
    b(1:N-2)=1+2*mu*TT;
    c(1:N-3)=-mu*TT;
    d(2:N-3)=transpose((1-2*mu*(1-TT))*U(3:N-2,k)+mu*(1-TT)*(U(2:N-3,k)+U(4:N-1,k)))+(1-TT)*r*f(X(3:N-2),T(1,k))+TT*r*f(X(3:N-2),T(1,k+1));
    U(2:N-1,k+1)=Thomas_algorithm(a,b,c,d);
end

figure
legend
for i=0:4
    plot(X(1:N),U(1:N,floor(M/5*i)+1),'DisplayName',strcat('T=',string(T(1,floor(M/5*i)+1))))
    hold on
end
plot(X(1:N),U(1:N,M-1),'DisplayName',strcat('T=',string(tspan(2))))
hold on
legend()
xlabel('X')
ylabel('U(T,X)')
hold off

figure
mesh(X,T,U);
xlabel('X')
ylabel('T')
zlabel('U')




