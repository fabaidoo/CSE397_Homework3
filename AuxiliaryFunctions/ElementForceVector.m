function Fe = ElementForceVector(quad_rule, P_e_mat, C_e, w_e, span1, span2)

%quadrature points and weights
[weights, xipoints] = gauss_legendre1D(quad_rule);

%pressure
pressure = 1;

Fe_partial = zeros(2*length(w_e), length(weights));

for i = 1:length(weights)
    xi1 = xipoints(i, :);
    weight = weights(i);
    
    %Jacobians
    [~, ~, ~, dxdv_phi] ...
        = ElementShapeFunctionRoutine([xi1 -1], P_e_mat, C_e, w_e, span1, span2);
    J_Gamma = norm(dxdv_phi(:, 1));
    
    JBez = Parent2BezierMap(span1, span2);
    JBez_gamma = JBez(1, 1);
    
    %Outward pointing normal
    nhat = - dxdv_phi(:, 2) / norm(dxdv_phi(:, 2));
    
    NMat = NMatrix([xi1 -1] , C_e, w_e, span1, span2);
    
    Fe_partial(:, i) = ...
        - weight * pressure * NMat' * nhat * J_Gamma * JBez_gamma; 
end

Fe = sum(Fe_partial, 2);

%--------------------------------------------------------------------------
function [weights, points] = gauss_legendre1D(quad)
        %produces the 1D gauss points and weights for a given "quad" rule
        %default is 2point rule
        
        if quad ~= 2 && quad ~= 3
           warning(...
               'Only 2 or 3 point quadrature allowed. Defaulting to 2 point...') 
        end
        
        if quad == 3
            points =...
                sqrt(3/5) * [-1; 0; 1];
            weights = [5 8 5]/ 9;
        else
            points = 1/ sqrt(3) * [-1 ; 1];
            weights = ones(1, 2);
        end
end

end