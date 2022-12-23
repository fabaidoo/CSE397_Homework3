function [DR_e_phi, R_e_phi] = ...
    ElementNURBSderivative(xi, C_e, w_e, span1, span2)
%Produces a matrix of the derivative of the NURBS basis at the bezier element 
%"span1 x span2" evaluated at "xi" in [-1, 1]^2. Uses the element 
%extraction operator "C_e" and the element NURBS weights "w_e" for this

size_we = size(w_e);
%Ensure vector of element weights is a column vector
if size_we(2) > size_we(1)
    w_e = w_e';  
end

p=2;
[Bplus, DBplus]= Bernstein2D(xi, p); %bernsteins on [-1, 1]^2
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

