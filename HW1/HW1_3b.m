
%% Problem 3 part (b)
clear all;clc;close all
% domain: x\in[a,b], y\in[c,d]
a = 0; b = 1;
c = 0; d = 1;
    
M_list = 10:10:100;
M_list = M_list+1;

for loop = 1:size(M_list,2)
    disp(loop)
    % grid sizes
    M = M_list(loop); 
    N = M; 
    hx = (b-a)/(M-1); 
    hy = (d-c)/(N-1); 
    h = hx;

    % 2D arrays of grids
    [X,Y] = meshgrid(a:hx:b,c-hy:hy:d+hy);

    % Initial f function
    f =  0.*X.*Y;

    % Initial setting: Boundary conditions
    u = 0.*X.*Y;
    u(:,1) = cos(2*pi*Y(:,1));
    u(1,1) = u(2,1);
    u(M+2,1) = u(M+1,1);
    uold = u;

    % Gauss-Seidel iteration 
    % reference code in Chapter 4.1 LeVeque, R. J. (2007).
    maxiter = 100000;
    
    for iter=0:maxiter
        unew = uold;
        unew(1,2:M-1) = unew(3,2:M-1);
        for j=2:(M-1)
            for i=2:(M+1)
                unew(i,j) = 0.25*(unew(i-1,j) + unew(i+1,j) + unew(i,j-1) + unew(i,j+1) - h^2 * f(i,j));
            end
        end
        unew(M+2,2:M-1) = unew(M,2:M-1);
        r = norm(unew-uold);
        if r < 10^(-6)
            break
        end
        uold = unew;
    end
    u_list{loop} = unew;
end

save HW1_3b.mat M_list u_list
        
%% part (b) plot
% reference professor's plot script
% plot of the maximum norm of the error vs.h, the grid spacing
clear all;clc;close all

load HW1_3b.mat
h = 1./(M_list-1);
error_list = zeros(size(M_list));

% calculate the maximum norm
for i = 1:size(u_list,2)
    exact_u = exact_sol(M_list(i));
    error = u_list{i} - exact_u;
    error_list(i) = norm(error,'Inf');
end

figure(1); clf();

loglog(h,error_list,'o-', 'LineWidth', 2)
hold on; 
loglog(h, h.^2, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell_\infty$ error','Interpreter','latex', 'FontSize', 24)

lgd = legend("error", "$\mathcal{O}(h^2)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';

saveas(gcf,'HW1_3b.png')
