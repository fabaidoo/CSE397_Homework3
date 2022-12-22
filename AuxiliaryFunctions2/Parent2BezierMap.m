function JMatrix = Parent2BezierMap(span1, span2)
% Returns the constant jacobian matrix of the affine map from the parent 
% element [-1, 1]^2 to the bezier element [a1, b1]x[a2, b2]. 

assert(length(span1) == 2 && length(span2) == 2,...
    'each span must be a 1D interval')
assert( span1(1) < span1(2) && span2(1) < span2(2),...
    'incorrectly defined interval')

a1 = span1(1);
b1 = span1(2);
a2 = span2(1);
b2 = span2(2);

J1 = 0.5 * ( b1 - a1);
J2 = 0.5 * (b2 - a2);

JMatrix = diag([J1 J2]);
end