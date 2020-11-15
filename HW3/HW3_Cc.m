% reference Trefethen Book p20.m - 2nd-order wave eq. on Chebyshev grid (compare p6.m)

% Time-stepping by leap frog formula:
clear all;clc;close all
N_list = 2.^[3:1:7];
for i = 1:length(N_list)
    N = N_list(i);
    disp(i)
    x = cos(pi*(0:N)/(N)); y = x';
    dt = 8/N^2;
    [xx,yy] = meshgrid(x,y);
    tmax = round(0.7/dt);

    % ut0
    vvt =exp(-100*(xx).^2).*exp(-100*(yy).^2); 
    %u0
    vvold = zeros(size(xx)); 
    %u1
    vv = vvold + dt*vvt  - 1/2*dt^2 * delt_op(vvold,x,y)+1/6*dt^3*delt_op(vvt,x,y)-1/24*dt^4*delt_op(delt_op(vvold,x,y),x,y)+dt^2*(delt_op(vvold,x,y)+1/12*dt^2*delt_op(delt_op(vvold,x,y),x,y));

    for j = 1:tmax
        vvnew = 2*vv - vvold + dt^2*(delt_op(vv,x,y)+1/12*dt^2*delt_op(delt_op(vv,x,y),x,y));
        
        vvold = vv; vv = vvnew;
        
        % view evolution
        figure(1);surf(xx,yy,vv);axis([-1 1  -1 1  -0.05 0.05]);shading interp;drawnow
    end
    u_list{i} = vv;
end

%% calculate error
for i = 1:length(N_list)-1
    N =  N_list(i);
    x = cos(pi*(0:N)/(N)); y = x';
    [xx,yy] = meshgrid(x,y);
    N_exact = N_list(end);
    x_exact = cos(pi*(0:N_exact)/(N_exact)); y_exact = x_exact';
    [xx_exact,yy_exact] = meshgrid(x_exact,y_exact);

    interp_exactU = interp2(xx_exact,yy_exact,u_list{end},xx,yy);
    error_list{i} = norm(u_list{i}-interp_exactU,'inf');
end
save HW3_Cc.mat N_list error_list
% plot results
clear all;clc;close all

load HW3_Cc.mat
N_list(end)=[];
h = 1./N_list;
a = [error_list{:}];

figure(1); clf();

loglog(h,a,'o-', 'LineWidth', 2)
hold on; 

loglog(h, h.^4, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell_\infty$ error','Interpreter','latex', 'FontSize', 24)

lgd = legend("error", "$\mathcal{O}(h^4)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';

saveas(gcf,'HW3_Cc.png') 
  
%% function
% CHEBFFT  Chebyshev differentiation via FFT. Simple, not optimal.  
%          If v is complex, delete "real" commands.

  function delta = delt_op(vv,x,y)
N =size(vv,1)-1;
uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
ii = 2:N;
for i = 2:N                % 2nd derivs wrt x in each row
  v = vv(i,:); V = [v fliplr(v(ii))];
  U = real(fft(V));
  W1 = real(ifft(1i*[0:N-1 0 1-N:-1].*U)); % diff wrt theta
  W2 = real(ifft(-[0:N 1-N:-1].^2.*U));    % diff^2 wrt theta
  uxx(i,ii) = W2(ii)./(1-x(ii).^2) - x(ii).* ... 
                 W1(ii)./(1-x(ii).^2).^(3/2);
end
for j = 2:N                % 2nd derivs wrt y in each column
  v = vv(:,j); V = [v; flipud(v(ii))];
  U = real(fft(V));
  W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));% diff wrt theta   
  W2 = real(ifft(-[0:N 1-N:-1]'.^2.*U));   % diff^2 wrt theta
  uyy(ii,j) = W2(ii)./(1-y(ii).^2) - y(ii).* ...
                 W1(ii)./(1-y(ii).^2).^(3/2);
end
delta = uxx+uyy;
  end