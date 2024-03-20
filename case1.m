% Initializing variables and parameters
dx=.05;
sr=1;
% dx=2^-9;
x=-5:dx:5; 
c=1;
dt=dx*sr/c; 
t= 0:dt:20;
u=zeros(length(t),length(x));

% Setting boundary conditions and initial conditions
% Dirichlet BCs sets first and last row of U to zero
% I.C f(x) = 1-|x| for |x| <=1

u(:,1)=0;
u(:,end)=0;
lx=length(x);
for i=1:lx
    if abs(x(i))<=1
        u(1,i)=1-abs(x(i));
    end
end


% u=zeros(length(t)+1,length(x));
s=c*dt/dx;
u(2,2:lx-1)=s^2/2*u(1,3:lx)+u(1,2:lx-1)*(1-s^2)+s^2/2*u(1,1:lx-2);
% u(1,:)=u(3,:);
% Solving for solution at other times
for i=3:length(t)
    u(i,2:lx-1)=s^2*u(i-1,3:lx)+u((i-1),(2:lx-1))*(2-2*s^2)+s^2*u((i-1),(1:lx-2))-u((i-2),(2:lx-1));
end
%%% plotting solution at different times
figure;plot(x,u(1,:))
t1=round(t,3);
hold on
plot(x,u(find(t1==5),:))
plot(x,u(find(t1==20),:)); 
hold off;
figure; plot(t,u(1:end,find(x==4.9)))
figure;surf(x,t1,u(1:end,:));colormap jet;shading interp