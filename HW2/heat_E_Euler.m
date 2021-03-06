function heat_E_Euler(a)
% 1D heat equation with Dirichlet boundary conditions
% u_t = u_xx
% Euler explicit ("Simplest")

N = 2^8;
h = 1/N;
x = (0:h:1)';
u0 = 10 * (1-x) .* exp(-200*(.5-x).^3) .* sin(30.*exp(-100*(x-.2).^2));
figure(1); clf; plot(x,u0);
pause

% final time
T = 0.01;

% time step
dt = a*h^2/2;       % choose a good time step (for stability)
r = dt/(2*h^2);
% indixes
I = 2:N;
% initial condition
U = u0;

% number of time steps
n_it = ceil(T/dt);
fprintf("Total number of time steps %d \n", n_it)

for n = 1:n_it
    
    U(I) = U(I) + 2*r*(U(I-1) - 2*U(I) + U(I+1));
    
    figure(1); clf;
    plot(x,U); axis([0 1 -3 5]);
    
end

end