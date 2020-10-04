%% Problem 3 part (f) 
% number of iteration for convergence

% domain: x\in[a,b], y\in[c,d]
a = 0; b = 1;
c = 0; d = 1;

% M_list and noiter_list
M_list = 10:10:100;
M_list = M_list+1;
noiter_list = zeros(size(M_list));

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
    
    % error check
    exact_u = exact_sol(M);
    e0 = norm(u - exact_u);
    epsilon = 50 * h^2;
    ek = e0 * epsilon;
    error = e0;
    
    % Gauss-Seidel iteration 
    % reference code in Chapter 4.1 LeVeque, R. J. (2007).
    noiter = 0;
    while error > ek
        noiter = noiter + 1;
        u(1,2:M-1) = u(3,2:M-1);
        for j=2:(M-1)
            for i=2:(M+1)
                u(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - h^2 * f(i,j));
            end
        end
        u(M+2,2:M-1) = u(M,2:M-1);
        current_error = u-exact_u;
        error = norm(u-exact_u);
    end
    noiter_list(loop) = noiter;
end

save HW1_3f.mat M_list noiter_list
%% plot the result
clear all;clc

load HW1_3f.mat

figure(1); clf();

hold on
plot(M_list.^2.*log(M_list),noiter_list,'LineWidth',2);
scatter(M_list.^2.*log(M_list),noiter_list,150,'filled');
hold off

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('$m^2logm$','Interpreter','latex', 'FontSize', 24)
ylabel('Number of iterations','Interpreter','latex', 'FontSize', 24)

saveas(gcf,'HW1_3f.png')


