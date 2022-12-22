function  [Bplus, DBplus]= Bernstein2D2(xi, p)
%Outputs a (p+1)x1 vector of the 2D degree p Bernstein basis at a point xi 
%on the standard element [-1, 1]^2 along with a (p+1)x2 matrix with their 
%gradient

%kron reversed


assert(length(xi)== 2 , 'Can only evaluate 2D points')
assert(xi(1) <= 1 && xi(1) >= -1 && xi(2) <= 1 && xi(2) >= -1,...
    'Evaluation point not in parent element [-1, 1]^2' )

xi1 = xi(1);
xi2 = xi(2);

%2D bernsteins go here 
%Bplus = zeros((p+1)^2, 1); %(Column vector!!) 
DBplus = zeros((p+1)^2, 2);

%1D bernsteins go here
B1 = zeros(p+1, 1); 
dB1 = zeros(p+1, 1);
B2 = B1;
dB2 = dB1;
for i = 1 : p+1
    %fill up 1D bernsteins
    [B1(i), dB1(i)] = Bernstein1D(xi1, i, p);
    [B2(i), dB2(i)] = Bernstein1D(xi2, i, p);
end

%Fill up the matrices
Bplus = kron(B2, B1);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DBplus(:, 1) = kron(B2, dB1); 
DBplus(:, 2) = kron(dB2, B1);


%------------------%------------------%------------------%-----------------
%Supplementary 1D Bernstein Function

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

%plot of bases for sanity check
if Flag == true 
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

end


