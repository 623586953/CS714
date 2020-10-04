
%% Problem 3 part (g)
clear all;clc;close all
% domain: x\in[a,b], y\in[c,d]
a = 0; b = 1;
c = 0; d = 1;

M_list = 10:10:100+1;
iteration_list = zeros(2,3);

for func = 1:2
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
        if func == 1
            u(:,1) = cos(2*pi*Y(:,1));
        else
            u(:,1) = sign(cos(2*pi*Y(:,1)));
        end
        u(1,1) = u(2,1);
        u(M+2,1) = u(M+1,1);
        uold = u;

        % Gauss-Seidel iteration 
        % reference code in Chapter 4.1 LeVeque, R. J. (2007).
        maxiter = 10000*M;
        r=[];
        no_iteration = 0;
        for iter=0:maxiter
            no_iteration = no_iteration +1;
            
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
            
            % update uold
            uold = unew;
        end
        iteration_list(func,loop) = no_iteration;
    end
end

save HW1_3g.mat M_list iteration_list
 %% plot the result
clear all;clc

load HW1_3g.mat

figure(1); clf();

bar(M_list,iteration_list.');

legend('','','Location','northwest');
set(gca,'FontSize',30,'box','on');

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('$M$','Interpreter','latex', 'FontSize', 24)
ylabel('Number of iterations','Interpreter','latex', 'FontSize', 24)
lgd = legend("$f(y) = cos(2\pi y)$", "$\hat{f}(y) = sgn\left(cos(2\pi y)\right)$",'FontSize', 24,...
       'Interpreter','latex');
   
saveas(gcf,'HW1_3g.png')