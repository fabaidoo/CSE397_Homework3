function dRe_dx= PhysicalDerivative(xi, P_e_mat, Ce, we, span1, span2)
% Gives the derivative of the element NURBS basis with respect to the
% physical coordinates

[inv_matrix, ~ , ~, ~] = ...
     ElementShapeFunctionRoutine(xi, P_e_mat, Ce, we, span1, span2);

[DR_e_phi, ~] = ...
    ElementNURBSderivative(xi, Ce, we, span1, span2);

dRe_dx = DR_e_phi * inv_matrix;

end