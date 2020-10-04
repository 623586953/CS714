% exact solution for problem c in HW1
function exact_u = exact_sol(M)
% domain: x\in[a,b], y\in[c,d]
a = 0; b = 1;
c = 0; d = 1;

N = M; 
hx = (b-a)/(M-1); 
hy = (d-c)/(N-1); 
h = hx;

% 2D arrays of grids
[X,Y] = meshgrid(a:hx:b,c-hy:hy:d+hy);
exact_u = cos(2*pi*Y).*(-1/(exp(4*pi)-1)*exp(2*pi*X)+exp(4*pi)/(exp(4*pi)-1)*exp(-2*pi*X));
end