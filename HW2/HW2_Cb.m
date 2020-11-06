
% HW 2 problem C.a
% reference Professor LEONARDO's code
% 2D heat equation with Dirichlet boundary conditions
% u_tt = u_xx+u_yy
% Euler explicit ("Simplest")
a = 0.9;
N_list = [2^10,2^9,2^8];
error_list = zeros(length(N_list)-1,1);
for i = 1:length(N_list)
    h = 1/N;
    [X,Y] = meshgrid(0:h:1,0:h:1);
    u0 = X.*Y.*0;

    % final time
    T = 0.01;

    % time step
    dt = a*h^2/2;       % choose a good time step (for stability)
    r = dt^2/(h^2);
    % indixes
    I = 2:N;
    % initial condition
    U = u0;
    Uold = -dt .* exp(-400*(X-0.5).^2).* exp(-400*(Y-0.5).^2);

    % number of time steps
    n_it = ceil(T/dt);
    fprintf("Total number of time steps %d \n", n_it)

    for n = 1:n_it
        disp(n/n_it)
        temp = U;
        U(I,I) = 2*U(I,I) - Uold(I,I) + r*(U(I-1,I)+U(I+1,I)+U(I,I-1)+U(I,I+1) - 4*U(I,I));
        Uold = temp;
    end
    
    if i == 1
        exact_sol = U;
    else
        error_list(i-1) = norm(U-exact_sol,inf);
    end
end
save HW2_Ca.mat N_list error_list
%% plot results
clear all;clc;close all

load HW2_Ca.mat
N_list(1) = [];
h = 1./N_list;
error_list

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

saveas(gcf,'HW3_.png')
