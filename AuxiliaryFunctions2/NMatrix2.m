function NMat =...
    NMatrix2(xi, C_e, w_e, span1, span2)
%Produces the strain rate tensor basis matrix evaluated at xi \in [-1, 1]

[~, R_e_phi] = ...
    ElementNURBSderivative2(xi, C_e, w_e, span1, span2);

NMat = zeros(2, 2*length(w_e));

for i = 1: length(w_e)
    Rei = R_e_phi(i);
    NMat(:, 2*(i-1)+1:2*i) = [Rei 0; 0 Rei];   
end

end