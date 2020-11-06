% HW2.b
clear all;clc;
N_list =  100:1:200;

for i = 1:length(N_list)
    %exact_folution
    N = N_list(i);
    x = 0:1/N:1;
    fx = exp(-400*(x-0.5).^2);
    X = 0:1/(N+1):1;
    fX = exp(-400*(X-0.5).^2);
    interp_fX = interp1(x,fx,X);
    error = norm(fX - interp_fX,inf);
    if abs(error) < 10^(-2)
        break
    end
end
disp(N_list(i))
%%

xx = 0.01;
100*(2*xx^3+xx^2)/(xx+1)^2