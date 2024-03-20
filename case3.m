% Initializing variables and parameters
dx=0.05;
% dx=2^-9;
x=-5:dx:5; 
c=1;
dt=dx/c; t= 0:dt:20;
u=zeros(length(t),length(x)+1);
% Setting boundary conditions and initial conditions
% Dirichlet BCs sets first and last row of U to zero
% I.C f(x) = 1-|x| for |x| <=1
u(:,1)=0;

lx=length(x);

for i=1:lx
    if abs(x(i))<=1
        u(1,i)=1-abs(x(i));
    end
end
u(1,end)=u(1,end-2);

% figure; plot(x,u(2,1:end-1))
s=c*dt/dx;
u(2,2:lx)=s^2/2*u(1,3:lx+1)+u(1,2:lx)*(1-s^2)+s^2/2*u(1,1:lx-1);  
% u(1,:)=u(3,:);
% Solving for solution at other times
for i=3:length(t)
 
    u(i,2:lx-1)=s^2*u(i-1,3:lx)+u((i-1),(2:lx-1))*(2-2*s^2)+s^2*u((i-1),(1:lx-2))-u((i-2),(2:lx-1));
    u(i-1,end)=-dx/dt*(u(i,end-2)-u(i-2,end-2))+u(i-1,end-2); 
    u(i,end-1)=s^2*u(i-1,end)+u(i-1,end-1)*(2-2*s^2)+s^2*u(i-1,end-2)-u(i-2,end-1);
    
%     u(i,2:lx)=s^2*u(i-1,3:lx+1)+u((i-1),(2:lx))*(2-2*s^2)+s^2*u((i-1),(1:lx-1))-u((i-2),(2:lx));
%     u(i,end)=u(i,end-2);
end
U=u(:,1:end-1);
%%% plotting solution at different times
figure;plot(x,U(1,:))
t1=round(t,3);
% hold on
% plot(x,U(find(t1==3),:))
% plot(x,U(find(t1==20),:)); 
% hold off;
% figure; plot(t,U(2:end,find(x==5)))
figure;surf(x,t,U(1:end,:));colormap jet;shading interp