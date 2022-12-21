function [e1,e2] = ElementNumber1D(e, n_el1)
%Given an element number, the functions below return the
%corresponding 1D element numbers of the 1D Bezier spans that the 2D element
%is made of. We assume that the formula used for the 2D element number is
%that the formula used for the global node number is e = n_el1(e2-1) + e1.

e1_prelim = mod(e, n_el1); %gives correct e1 except when e= m*n_el1 then equals 0
if e1_prelim == 0
    e1 = n_el1;
else
    e1 = e1_prelim;
end
e2 = (e - e1)/ n_el1 + 1;
end

