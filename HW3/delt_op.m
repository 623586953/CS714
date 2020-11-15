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