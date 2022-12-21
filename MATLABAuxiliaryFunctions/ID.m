function [P, ArrayMatrix] = ID(l, B)
%Credit to Minyeong Kim for this code
%ID array
% input  
%   i: degree of freedom
%   A: global node number
% output
%  P: Global equation number

n_1 = 5;
n_2 = 5;
IDArray = zeros(2, n_1*n_2);

k = 1;
for j = 1:n_2
    for i = 1:n_1
        A = n_1*(j-1)+i;
        if i == 1
            IDArray(1, A) = 0;
            IDArray(2, A) = k;
            k = k+1;
        elseif i == n_1
            IDArray(1, A) = k;
            IDArray(2, A) = 0;
            k = k+1;
        else
            IDArray(1, A) = k;
            IDArray(2, A) = k+1;
            k = k+2;
        end
    end
end

P = IDArray(l, B);


    ArrayMatrix = zeros(2, n_1*n_2);
    for n=1:2
        for G = 1: n_1*n_2
           ArrayMatrix(n, G)  = IDArray(n, G); 
        end  
    end
   

end

