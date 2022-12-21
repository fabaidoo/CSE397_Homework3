function Ke = ElementStiffness(quad_rule, P_e_mat, C_e, w_e, span1, span2)
%Builds the element stiffness matrix based of a 2D gauss quadrature
%determined by quad. Only 2point and 3 point quadrature is possible.

%quadrature points and weights
[weights, xipoints] = gauss_legendre2D(quad_rule);

%D matrix
YoungModulus = 1;
PoissonRatio = 0.3;
DMat = YoungModulus /(1 - PoissonRatio^2) * ...
    [1 PoissonRatio 0; PoissonRatio 1 0 ; 0 0 (1 - PoissonRatio)/2];

%contribution of each quad point will go here 
Ke_partial = zeros(2*length(w_e), 2*length(w_e), length(weights)); 

for i = 1:length(weights)
    xi = xipoints(i, :);
    weight = weights(i);
    
    %Jacobians
    [~, Jx, ~, ~] ...
        = ElementShapeFunctionRoutine(xi, P_e_mat, C_e, w_e, span1, span2);
    Jphi = Parent2BezierMap(span1, span2);
    
    %BMatrix
    BMat =...
    BMatrix(xi, P_e_mat, C_e, w_e, span1, span2);

    %Quad point contribution to element stiffness matrix
    Ke_partial(:,:, i) = weight * BMat' * DMat * BMat * det(Jphi) * abs(Jx);    
end

%Add up contributions to get the element stiffness matrix
Ke = sum(Ke_partial, 3);


%-------------------%----------------------%-------------------------%-----
    function [weights, points] = gauss_legendre2D(quad)
        %produces the 2D gauss points and weights for a given "quad" rule
        %default is 2point rule
        
        if quad ~= 2 && quad ~= 3
           warning(...
               'Only 2 or 3 point quadrature allowed. Defaulting to 2 point...') 
        end
            
        if quad == 3
            points =...
                sqrt(3/5) * [0 0; 0 -1; 0 1; -1 0;...
                -1 -1; -1 1; 1 0; 1 -1; 1 1];
            weights = [64 40 40 40 25 25 40 25 25]/ 81;
        else
            points = 1/ sqrt(3) * [-1 -1; -1 1; 1 -1; 1 1];
            weights = ones(1, 4);
        end
 
    end




end