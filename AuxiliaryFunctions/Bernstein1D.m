function [B, dB] = Bernstein1D(xi, i, p, Flag)
% The i-th degree p 1D Bernstein basis function evaluated at xi in [-1, 1]
% along with its derivative
assert(p >= 0, 'Polynomial Degree must be non-negative')
assert( xi >= -1 && xi <= 1, 'Can only evaluate within [-1, 1]' )
if nargin == 3
    Flag = false ;
end

%recursive definition of bernstein polynomials
if p == 0
    if i < 1 || i > p+1
        B = 0;
        dB = 0;
    else
        B = 1;
        dB = 0;
    end
    
else
    [Btemp1, ~]= Bernstein1D(xi, i, p-1);
    [Btemp2, ~] = Bernstein1D(xi, i-1, p-1);
    B = 0.5 * ((1 - xi)*Btemp1 + (1 + xi) * Btemp2) ;
    dB = 0.5 * p * (Btemp2 - Btemp1);   
end

if Flag == true %plot of basis
   x = linspace(-1, 1, 100);
   Btest = zeros(length(x), p+1);
   dBtest = zeros(length(x), p+1);
   fig1 = figure;
   fig2 = figure;
   for k = 1: p+1
       for l = 1: length(x)
            [Btest(l, k), dBtest(l, k)] = Bernstein1D(x(l), k, p);
       end
   
    figure(fig1)
    plot(x', Btest(:, k) , 'LineWidth', 1)
    hold on
    
    figure(fig2)
    plot(x', dBtest(:, k), '--', 'LineWidth', 1)
    hold on
   end     
end
end