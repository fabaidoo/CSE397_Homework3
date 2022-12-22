function [inv_matrix, Jacobian, x_phi, dxdv_phi] ...
    = ElementShapeFunctionRoutine2(xi, P_e_mat, Ce, we, span1, span2)
%Element shape function routine outputs the inverse of the Jacobian matrix
%and the (absolute value) of the Jacobian determinant evaluated at a point
%xi in [-1, 1]^2

sizechk = size(P_e_mat);
assert(sizechk(2) == 2, 'Matrix of control points should be n_e x 2')

[DRephi, Rephi] = ElementNURBSderivative2(xi, Ce, we, span1, span2); 
x_phi = P_e_mat' * Rephi;
dxdv_phi = P_e_mat' * DRephi;

inv_matrix = inv(dxdv_phi);
Jacobian = (det(dxdv_phi));


function [DR_e_phi, R_e_phi] = ...
    ElementNURBSderivative2(xii, C_e, w_e, span1, span2)
%Produces a matrix of the derivative of the NURBS basis at the bezier element 
%"span1 x span2" evaluated at "xi" in [-1, 1]^2. Uses the element 
%extraction operator "C_e" and the element NURBS weights "w_e" for this

size_we = size(w_e);
%Ensure vector of element weights is a column vector
if size_we(2) > size_we(1)
    w_e = w_e';  
end

p=2;
[Bplus, DBplus]= Bernstein2D2(xii, p); %bernsteins on [-1, 1]^2
Jphi = Parent2BezierMap(span1, span2); %jacobian matrix of phi^e map

W_e = diag(w_e); %diagonal matrix of weights
%weighting function and its derivative
wv = w_e' * (C_e * Bplus);
dwv = w_e' * (C_e * DBplus/Jphi );
%dwv = w_e' * (C_e * DBplus * inv(Jphi) );

%vector of NURBS
R_e_phi = W_e * C_e * Bplus / wv ;
%matrix of NURBS derivatives
DR_e_phi = W_e * C_e * (DBplus/Jphi - Bplus * dwv/wv ) / wv ;
end

end