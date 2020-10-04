
%% HW1-D
% regular 
clear all;clc;close all
% domain: x\in[a,b], y\in[c,d]
a = 0; b = 1;
c = 0; d = 1;

% grid sizes
M = 101;
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
maxiter = 100000*M;

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

exact_u = exact_sol(M);
error = norm(exact_u - unew);
% Show results
disp(sprintf('Take %d steps to reach convergence, and the error is %d',iter,error))

        
%% Multigrid 
% h to 2h to h
clear all;clc;close all
% domain: x\in[a,b], y\in[c,d]
a = 0; b = 1;
c = 0; d = 1;

% grid sizes
M = 101;
N = M; 
hx = (b-a)/(M-1); 
hy = (d-c)/(N-1); 
h = hx;

% 2D arrays of grids
[X,Y] = meshgrid(a:hx:b,c-hy:hy:d+hy);
[X2h,Y2h] = meshgrid(a:hx*2:b,c-hy*2:hy*2:d+hy*2);

% Initial f function
f =  0.*X.*Y;

% Initial setting: Boundary conditions
u = 0.*X.*Y;
u(:,1) = cos(2*pi*Y(:,1));
u(1,1) = u(2,1);
u(M+2,1) = u(M+1,1);
uold = u;

% 1. Take 3 iterations on initial grid.
for iter = 1:3
    unew = uold;
    unew(1,2:M-1) = unew(3,2:M-1);
    for j=2:(M-1)
        for i=2:(M+1)
            unew(i,j) = 0.25*(unew(i-1,j) + unew(i+1,j) + unew(i,j-1) + unew(i,j+1) - h^2 * f(i,j));
        end
    end
    unew(M+2,2:M-1) = unew(M,2:M-1);
    uold = unew;
end

% 2. Compute the residual, and then coarsen the residual.\\
residual = zeros(size(unew));
for j=2:(M-1)
    for i=2:(M+1)
        residual(i,j) = (4* unew(i,j) - unew(i-1,j) - unew(i+1,j)- unew(i,j-1) - unew(i,j+1))/h^2 + f(i,j);
    end
end

residual_2h = residual;
residual_2h(3:2:end-2,:) = [];
residual_2h(:,2:2:end) = [];
residual_2h(1,:) = residual_2h(3,:);
residual_2h(size(residual_2h,1),:) = residual_2h(size(residual_2h,1)-2,:);

% 3. Approximate the error on the coarse grid, and fit it back to the initial grid.\\
enew_2h = zeros(size(residual_2h));
re = 1;

while re > 0.04
    eold_2h = enew_2h;
    enew_2h(1,2:(M-1)/2) = enew_2h(3,2:(M-1)/2);
    for j=2:size(enew_2h,2)-1
        for i=2:size(enew_2h,1)-1
            enew_2h(i,j) = 0.25*(enew_2h(i-1,j) + enew_2h(i+1,j) + enew_2h(i,j-1) + enew_2h(i,j+1) - (2*h)^2 * residual_2h(i,j));
        end
    end
    enew_2h(size(enew_2h,1),2:size(enew_2h,2)-1) = enew_2h(size(enew_2h,1)-2,2:size(enew_2h,2)-1);
    re = norm(enew_2h- eold_2h);
end

enew_h = interp2(enew_2h,1);
enew_h(size(enew_h,1),:) = [];
enew_h(1,:) = [];

% 4. Update the approximation on the initial grid, and iterate until converge.
unew = unew - enew_h;
uold = unew;

for iter2 = 1:10^6
    unew = uold;
    unew(1,2:M-1) = unew(3,2:M-1);
    for j=2:(M-1)
        for i=2:(M+1)
            unew(i,j) = 0.25*(unew(i-1,j) + unew(i+1,j) + unew(i,j-1) + unew(i,j+1) - h^2 * f(i,j));
        end
    end
    unew(M+2,2:M-1) = unew(M,2:M-1);
    r = norm(unew-uold);
    uold = unew;
    if r < 10^(-6)
        break
    end
end

exact_u = exact_sol(M);
error =norm(exact_u-unew);

% Show results
disp(sprintf('Take %d steps to reach convergence, and the error is %d',iter2+3,error))
