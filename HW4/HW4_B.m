%% HW4B Godunov shceme

%% singel case -- using the conservative form of the equation

clear all;clc;close all;
    
% Grid and initial data:
N = 2^8;
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
    unew(2:N) = uold(2:N) - dt/h*(nf(uold(2:N),uold(3:N+1)) - nf(uold(1:N-1),uold(2:N)));
    unew(1) = uold(1) - dt/h*(nf(uold(1),uold(2))-nf(uold(N+1),uold(1)));
    unew(N+1) = uold(N+1) - dt/h*(nf(uold(N+1),uold(1))-nf(uold(N),uold(N+1)));
     % view evolution
    figure(1);plot(x,u0,x,unew);axis([0 1  0.5 2.5]);shading interp;drawnow
end

%% function
%numerical flux function: nf
function result = nf( u, v )
  for i = 1:length(u)
    if (u(i) >= v(i))
      % direction of shock speed
      if ((u(i)+v(i))/2 > 0)
        scriptu(i)=u(i);
      else
        scriptu(i)=v(i);
      end

    else

      if (u(i)>0)
        scriptu(i)=u(i);
      elseif (v(i)<0)
        scriptu(i)=v(i);
      else
        %stagnation condition
        scriptu(i)=0;
      end
    end
  end
  result = 0.5 * scriptu.^2;
  return
end
