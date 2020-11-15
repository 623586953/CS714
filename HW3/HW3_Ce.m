% Spectral method
clear all;clc;close all

B_list = 1:1:20;
parfor i = 1:length(B_list)
    disp(i)
    B = B_list(i);
    
    % Spectral method vv1
    N = 128;
    x = cos(pi*(0:N)/(N)); y = x';
    dt = 8/N^2;
    [xx,yy] = meshgrid(x,y);
    T = 0.75;
    maxiter = round(T/dt);

    % ut0
    vvt =sin(B*pi*(xx+1)/2).*sin(B*pi*(yy+1)/2); 
    % u0
    vvold = zeros(size(xx)); 
    % u1
    vv = vvold + dt*vvt  - 1/2*dt^2 * delt_op(vvold,x,y)+1/6*dt^3*delt_op(vvt,x,y)-1/24*dt^4*delt_op(delt_op(vvold,x,y),x,y)+dt^2*(delt_op(vvold,x,y)+1/12*dt^2*delt_op(delt_op(vvold,x,y),x,y));

    for j = 1:maxiter
        vvnew = 2*vv - vvold + dt^2*(delt_op(vv,x,y)+1/12*dt^2*delt_op(delt_op(vv,x,y),x,y));
        vvold = vv; vv = vvnew;
%         %plot3(xx,yy,vv)
%         % view evolution
%         if mod(j,10) == 0
%             figure(1);surf(xx,yy,vv);axis([-1 1  -1 1  -0.5 0.5]); shading interp;drawnow
%         end
    end
    xx1 = xx; yy1 = yy; vv1 = vv;
    
     % Spectral method vv2
    N = 128*2;
    x = cos(pi*(0:N)/(N)); y = x';
    dt = 8/N^2;
    [xx,yy] = meshgrid(x,y);
    T = 0.75;
    maxiter = round(T/dt);

    % ut0
    vvt =sin(B*pi*(xx+1)/2).*sin(B*pi*(yy+1)/2); 
    %u0
    vvold = zeros(size(xx)); 
    %u1
    vv = vvold + dt*vvt  - 1/2*dt^2 * delt_op(vvold,x,y)+1/6*dt^3*delt_op(vvt,x,y)-1/24*dt^4*delt_op(delt_op(vvold,x,y),x,y)+dt^2*(delt_op(vvold,x,y)+1/12*dt^2*delt_op(delt_op(vvold,x,y),x,y));

    for j = 1:maxiter
        vvnew = 2*vv - vvold + dt^2*(delt_op(vv,x,y)+1/12*dt^2*delt_op(delt_op(vv,x,y),x,y));
        vvold = vv; vv = vvnew;
%         %plot3(xx,yy,vv)
%         % view evolution
%         if mod(j,10) == 0
%             figure(1);surf(xx,yy,vv);axis([-1 1  -1 1  -0.5 0.5]); shading interp;drawnow
%         end
    end
    xx2 = xx; yy2 = yy; vv2 = vv;
    
    v_exact = interp2(xx2,yy2,vv2,xx1,yy1);
    SPerror_list(i) = norm(vv1-v_exact,'inf');
    
    N = 128;
    dt = 1/N^2;
    T = 0.75;
    maxiter = round(T/dt);
    
    %FD method
    u1 = wave_solution(B, N, T, dt, 0, 100);
    u2 = wave_solution(B, 2*N, T, dt/4, 0, 100);
    
    % grids in each dimension
    [x y] = meshgrid(-1:2/(N-1):1,-1:2/(N-1):1);
    [xx yy] = meshgrid(-1:2/(2*N-1):1,-1:2/(2*N-1):1);
    u_exact = interp2(xx,yy,u2,x,y);
    FDerror_list(i) = norm(u1-u_exact,'inf');
end

save HW3_ce.mat B_list FDerror_list SPerror_list

%%

%% plot results
clear all; clc;
load HW3_ce.mat
figure(1); clf();
hold on
semilogy(B_list,SPerror_list,'o-', 'LineWidth', 2);
semilogy(B_list,FDerror_list,'o-', 'LineWidth', 2);
plot([1 12],[0.03,0.03]);
hold off
ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell_\infty$ error','Interpreter','latex', 'FontSize', 24)
xlim([1,12]);
lgd = legend("Spectral Method", "Finite Difference Method",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';
%
saveas(gcf,'HW3_Ce.png') 

%% debug

% B_list = [1,2,3];
% for i = 1:length(B_list)
%     disp(i)
%     B = B_list(i);
%     N = 64;
%     dt = 1/N^2;
%     T = 1.;
%     maxiter = round(T/dt);
%     
%     %FD method
%     u1 = wave_solution(B, N, T, dt, 0, 100);
%     u2 = wave_solution(B, 2*N, T, dt/4, 0, 100);
%     
%     % grids in each dimension
%     [x y] = meshgrid(-1:2/(N-1):1,-1:2/(N-1):1);
%     [xx yy] = meshgrid(-1:2/(2*N-1):1,-1:2/(2*N-1):1);
%     u_exact = interp2(xx,yy,u2,x,y);
%     FDerror_list(i) = norm(u1-u_exact,'inf');
% end