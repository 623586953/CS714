%% sub-domain fixed length
% effect of overlapped region
clear all;clc;close all

index = 3;
a = 3;
b = 1;
n = 30*2^index;
m = 10*2^index;
dx = a/(n+1);
dy = b/(m+1);
x = 0:dx:a;
y = 0:dy:b;
[X Y] = meshgrid(x,y) ;

f1 = cos(2*pi*y);
f2 = cos(2*pi*y);

% initial conditon
domain_u0 = zeros(size(X));
domain_u0(:,1) = f1;
domain_u0(:,end) = f2;

overlap_list = 0.1:0.05:0.5;
iter_list = zeros(size(overlap_list,2),1);
for nnn = 1:size(overlap_list,2)
    % calculate the index of sub-domain
    num_subdomains = 3;
    overlap = overlap_list(nnn);
    subdomain_size = round((n+2)/(num_subdomains-overlap*(num_subdomains-1)));
    subdomain_index = zeros(num_subdomains,2);

    subdomain_index(1,1) = 1;
    subdomain_index(1,2) = subdomain_size;
    for i = 2:num_subdomains-1
        subdomain_index(i,1) = subdomain_index(i-1,2) - round(subdomain_size*overlap);
        subdomain_index(i,2) = subdomain_index(i,1) + subdomain_size;
    end
    subdomain_index(num_subdomains,1) = subdomain_index(num_subdomains-1,2) - round(subdomain_size*overlap);
    subdomain_index(num_subdomains,2) = (n+2);

    % Iteration
    domain_uold = domain_u0;
    for iter = 1:10^6
        % solve the problem in each subdomain
        for i = 1:num_subdomains
            subdomain_sol{i} = direct_solver(subdomain_index(i,2)-subdomain_index(i,1)-1,m,dx,dy,domain_uold(:,subdomain_index(i,1)),domain_uold(:,subdomain_index(i,2)));
        end
        % merge subdomain solution
        domain_unew = zeros(size(domain_uold));
        domain_unew(:,1) = f1;
        domain_unew(:,end) = f2;

        domain_unew(:,subdomain_index(1,1)+1:subdomain_index(1,2)-1) = subdomain_sol{1};
        for i = 2:num_subdomains
            temp = subdomain_sol{i};
            domain_unew(:,subdomain_index(i,1)+1:subdomain_index(i-1,2)-1) = (domain_unew(:,subdomain_index(i,1)+1:subdomain_index(i-1,2)-1) + temp(:,1:subdomain_index(i-1,2)-subdomain_index(i,1)-1))/2;
            domain_unew(:,subdomain_index(i-1,2):subdomain_index(i,2)-1)= temp(:,subdomain_index(i-1,2)-subdomain_index(i,1):end);
        end
        res = norm(domain_unew(:) - domain_uold(:))/norm(domain_uold(:));
        if res < 10^-6
            break
        end
        domain_uold = domain_unew;
    end
    iter_list(nnn) = iter;
end
    
%% plot result

figure(1); clf();
plot(overlap_list*2,iter_list,'o-', 'LineWidth', 2)

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('Overlapping area $(\%)$','Interpreter','latex', 'FontSize', 24)
ylabel('Iteration for convergence','Interpreter','latex', 'FontSize', 24)
    
saveas(figure(1),'Overlapping_area.png') 