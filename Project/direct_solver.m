% reference solution for HW1
% f1 is the left boundary condition
% f2 is the right boundary condition
function sol = direct_solver(n,m,dx,dy,f1,f2)

% Building the 1D Dirichelet minus Laplacian matrix
e = ones(n,1);
Asp = spdiags([-e 2*e -e], -1:1, n, n);

% Building the 1D Neumann minus Laplacian matrix
e = ones(m+2,1);
Bsp = spdiags([-e 2*e -e], -1:1, m+2, m+2);
Bsp(1,1:3) = [1.5 -2 0.5 ]*dy;
Bsp(end,end-2:end) = [ 0.5 -2 1.5 ]*dy;

% creating the identities (here carefull with the 
% boundaries)
I_A = speye(m+2,m+2);   I_A(1,1) = 0;   I_A(end,end) = 0;
I_B = speye(n,n);

% assembling the 2D minus Laplacian
Delta = kron(Asp/dx^2,I_A) + kron(I_B,Bsp/dy^2);

% writing the source term
%f1 = sin(2*pi*y);
f1(1) = 0;   f1(end) = 0;
e_1 = zeros(n,1);   e_1(1) = 1;
%f2 = sin(2*pi*y);
f2(1) = 0;   f2(end) = 0;
e_2 = zeros(n,1);   e_2(end)  = 1;

F = kron(e_1, f1).'+kron(e_2, f2).';
f = F(:)/dx^2;

% solving the system using a sparse solver for reference
u_dir = Delta\f;
sol = reshape(u_dir.',[],n);
end


