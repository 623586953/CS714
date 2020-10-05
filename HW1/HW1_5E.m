%% Problem 5E
% number of iteration for convergence
load HW1_5E_fine_grid.mat 
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
    u(:,1) = sign(cos(2*pi*Y(:,1)));
    u(1,1) = u(2,1);
    u(M+2,1) = u(M+1,1);
    
    % error check
    exact_u = interp2(finex,finey,finesol,X,Y);
    exact_u(1,:) = exact_u(3,:);
    exact_u(size(exact_u,1),:) = exact_u(size(exact_u,1)-2,:);
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
        if noiter > 100000
            break
        end
    end
    noiter_list(loop) = noiter;
end

save HW1_5E.mat M_list noiter_list
%% plot the result
clear all;clc


figure(1); clf();

hold on

load HW1_3f.mat
plot(M_list.^2.*log(M_list),noiter_list,'-o','LineWidth',2);
load HW1_5E.mat
plot(M_list.^2.*log(M_list),noiter_list,'-o','LineWidth',2);

hold off

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('$m^2logm$','Interpreter','latex', 'FontSize', 24)
ylabel('Number of iterations','Interpreter','latex', 'FontSize', 24)
lgd = legend("$f(y) = cos(2\pi y)$", "$\hat{f}(y) = sgn\left(cos(2\pi y)\right)$",'FontSize', 24,...
       'Interpreter','latex','Location','east');
   
saveas(gcf,'HW1_5E.png')


%% Problem 5E solution on find grid
% number of iteration for convergence

% domain: x\in[a,b], y\in[c,d]
a = 0; b = 1;
c = 0; d = 1;

% grid sizes
M = 400+1; 
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
u(:,1) = sign(cos(2*pi*Y(:,1)));
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

finesol = unew;
finex = X;
finey = Y;

save  HW1_5E_fine_grid.mat finesol finex finey


