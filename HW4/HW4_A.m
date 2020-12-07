%% HW4

%% singel case -- using the conservative form of the equation

clear all;clc;close all;
    
% Grid and initial data:
N = 2^7;
h = 1/N;
x = 0:h:1;

% final time
dt = h/4;
T = 0.7;

% number of iterations
Nit = round(T/dt);

% initial condition
u0 = 3/2+sin(2*pi*x);
unew = u0;

for i = 1:Nit
    uold = unew;
    unew(2:N) = uold(2:N) - 1/2*dt/h*(uold(2:N).^2-uold(1:N-1).^2);
    unew(1) = uold(1) - 1/2*dt/h*(uold(1).^2-uold(N+1).^2);
    unew(N+1) = uold(N+1) - 1/2*dt/h*(uold(N+1).^2-uold(N).^2);
    if mod(i,10) == 0
         % view evolution
        figure(1);plot(x,u0,x,unew);axis([0 1  0.5 2.5]);shading interp;drawnow
    end
end

%% singel case -- using the nonconservative form of the equation
clear all;clc;close all

% Grid and initial data:
N = 2^7;
h = 1/N;
x = 0:h:1;

% final time
dt = h/4;
T = 0.75;

% number of iterations
Nit = round(T/dt);

% initial condition
u0 = 3/2+sin(2*pi*x);
unew = u0;

for i = 1:Nit
    uold = unew;
    
    unew(N+1) = uold(N+1) - dt/h*((double(uold(N+1)>0).*uold(N+1)).*(uold(N+1)-uold(N))+(double(uold(N+1)<0).*uold(N+1)).*(uold(1)-uold(N+1)));
    unew(2:N) = uold(2:N) - dt/h*((double(uold(2:N)>0).*uold(2:N)).*(uold(2:N)-uold(1:N-1))+(double(uold(2:N)<0).*uold(2:N)).*(uold(3:N+1)-uold(2:N)));
    unew(1) = uold(1) - dt/h*((double(uold(1)>0).*uold(1)).*(uold(1)-uold(N+1))+(double(uold(1)<0).*uold(1)).*(uold(2)-uold(1)));
    if mod(i,10) == 0
         % view evolution
        figure(1);plot(x,u0,x,unew);axis([0 1  0.5 2.5]);shading interp;xlabel('$x$','Interpreter','latex', 'FontSize', 24);ylabel('u','Interpreter','latex', 'FontSize', 24);legend("Initial $u$", "numerical solution",'FontSize', 24,'Interpreter','latex');drawnow
    end
end

saveas(gcf,'HW4_numerical_solution.png') 

%% accuracy

% problem A
clear all;clc;close all;

N_list = 2.^[6:1:9];
for j = 1:length(N_list)
    disp(j)
    
    % Grid and initial data:
    N = N_list(j);
    h = 1/N;
    x = 0:h:1;


    % final time
    dt = h/4;
    T = 0.7;

    % number of iterations
    Nit = round(T/dt);

    % initial condition
    u0 = 3/2+sin(2*pi*x);
    unew = u0;

    for i = 1:Nit
        uold = unew;
        unew(N+1) = uold(N+1) - dt/h*((double(uold(N+1)>0).*uold(N+1)).*(uold(N+1)-uold(N))+(double(uold(N+1)<0).*uold(N+1)).*(uold(1)-uold(N+1)));
        unew(2:N) = uold(2:N) - dt/h*((double(uold(2:N)>0).*uold(2:N)).*(uold(2:N)-uold(1:N-1))+(double(uold(2:N)<0).*uold(2:N)).*(uold(3:N+1)-uold(2:N)));
        unew(1) = uold(1) - dt/h*((double(uold(1)>0).*uold(1)).*(uold(1)-uold(N+1))+(double(uold(1)<0).*uold(1)).*(uold(2)-uold(1)));
      end
    
    u_list{j} = unew;
end
%% plot accuracy
for i = 1:length(N_list)-1
    N =  N_list(i);
    x = 0:1/N:1; 
    N_exact = N_list(end);
    x_exact = 0:1/N_exact:1;

    interp_exactU = interp1(x_exact,u_list{end},x);
    error_list{i} = norm(u_list{i}-interp_exactU,'inf');
end

N_list(end)=[];
h = 1./N_list;
a = [error_list{:}];

figure(1); clf();

loglog(h,a,'o-', 'LineWidth', 2)
hold on; 

loglog(h, h, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell_\infty$ error','Interpreter','latex', 'FontSize', 24)

lgd = legend("error", "$\mathcal{O}(h)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'west';

saveas(gcf,'HW4_accuracy.png') 

%% Godunov shceme

% numerical flux function: nf