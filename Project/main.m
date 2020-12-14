%% Reference the solution of HW1

%% section 1. find the exact solution using exact solver
clear all;clc;close all

%initial settting
a = 3;
b = 1;
n = 300;
m = 100;
dx = a/(n+1);
dy = b/(m+1);
x = 0:dx:a;
y = 0:dy:b;
[X Y] = meshgrid(x,y) ;

f1 = cos(2*pi*y);
f2 = cos(2*pi*y);
%  director solver
sol_exact = direct_solver(300,100,dx,dy,f1,f2);

% % plot result
figure(1);surf(X(:,2:end-1),Y(:,2:end-1),sol_exact);xlabel('x');ylabel('y');shading interp;drawnow
figure(2);contourf(X(:,2:end-1),Y(:,2:end-1),sol_exact);colorbar('eastoutside');xlabel('x');ylabel('y');shading interp;drawnow
saveas(figure(1),'solution1.png') 
saveas(figure(2),'solution2.png') 

%% section 2
clear all;clc;close all
% number of refinements
Np = 5;
% arrays to store the error and mesh size
err = zeros(Np,1);
h = zeros(Np,1);

for index = 1:Np
    a = 1;
    b = 1;
    n = 10*2^index;
    m = 10*2^index;
    dx = a/(n+1);
    dy = b/(m+1);
    x = 0:dx:a;
    y = 0:dy:b;
    [X Y] = meshgrid(x,y) ;

    f1 = cos(2*pi*y);
    f2 = 0* y;
    
    sol = direct_solver(1,1,n,m,f1,f2);

    % sampling the analytical solution 
    sol_exact = cos(2*pi.*Y(:,2:end-1)).*(exp(-2*pi.*X(:,2:end-1))*exp(4*pi)- exp(2*pi.*X(:,2:end-1)))/( exp(4*pi)-1 );
    
    % computing the relative l^2 error
    err(index) = norm(sol_exact(:) - sol(:))/norm(sol_exact(:));
    % saving the step size
    h(index) = dx;
end

% Plotting the results 
figure(1); clf();
loglog(h,err,'o-', 'LineWidth', 2)
hold on; 
loglog(h, h.^2, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error', 'FontSize', 24);
xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^2$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h^2)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';
    
saveas(figure(1),'convergence_rate.png') 

%% section 3
% Effect of overlapped area on the convergence rate
a = 1;
b = 1;
n = 10*2^index;
m = 10*2^index;
dx = a/(n+1);
dy = b/(m+1);
x = 0:dx:a;
y = 0:dy:b;
[X Y] = meshgrid(x,y) ;

f1 = cos(2*pi*y);
f2 = 0* y;

% initial condition
u0 = zeros(size(X));
u0(:,1) = f1;
u0(:,end) = f2;

% calculate the domain of each sub-domain
num_subdomains = 3;
overlap = 0.2;


