function BMat =...
    BMatrix(xi, P_e_mat, C_e, w_e, span1, span2)
%Produces the strain rate tensor basis matrix evaluated at xi \in [-1, 1] 

dRe_dx = PhysicalDerivative(xi, P_e_mat, C_e, w_e, span1, span2);

BMat = zeros(3, 2*length(w_e));

for i = 1: length(w_e)
    dRei_dx = dRe_dx( i, :);
    BMat(:, 2*(i-1)+1:2*i) = [dRei_dx(1) 0; 0 dRei_dx(2) ; dRei_dx(2) dRei_dx(1)];   
end

end

