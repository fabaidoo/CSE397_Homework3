function [span1vec, span2vec] = BezierIntervals(Xi1, Xi2)
% Creates 1D bezier intervals in each direction using their knot vectors
%spanIvec is a n_elI x 2 matrix. Each row [a  b] represents a particular
%Bezier interval

%remove repeated knots
Xi1list = unique(Xi1);
Xi2list = unique(Xi2);

%preallocation
span1vec = zeros(length(Xi1list) - 1, 2);
span2vec = zeros(length(Xi2list) - 1, 2);

for i = 1:(length(Xi1list) - 1)
    span1vec(i, :) = Xi1list(i:i+1);   
end

%Because Xi1=Xi2 for our problem, I could have forego the loop below in favor of
%filling span2vec in the above for loop
for j = 1:(length(Xi2list) - 1)
    span2vec(j, :) = Xi2list(j:j+1);   
end

end