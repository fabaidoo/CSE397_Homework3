function displacements = implementation2(quad_rule)
%Solves for the displacements for homework 3. displacements have the same
%shape as the ID array. Specify the quadrature rule with the quad_rule
%input. Only 2 point and 3 point quadrature is possible. (You could add
%more quadrature rules in ElementForceVector.m AND ElementStiffness.m)


%IDArray
[~, IDArray] = ID(1,1);

%w: vector of weights
%P: matrix of control points. An (n_1*n_2) x 2 matrix
%Xi_1, Xi_2: Knot vectors
%p_1, p_2: basis degree
%n_1, n_2: number of bases in each direction
[w, Points, Xi_1, Xi_2, n_1, n_2, p_1, p_2] = DomainGeometry;

%n_el: number of elements
%C_operators: array of element extraction operators
% IEN: IEN array
[n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2);

%Matrix of bezier intervals in each direction
[span1vec, span2vec] = BezierIntervals(Xi_1, Xi_2);
n_el1 = length(span1vec);
%n_el2 = length(span2vec);

%quadrature rule
%quad_rule = 3;

%pre-allocation
neq = max(IDArray, [], 'all'); %number of equations
K = zeros(neq, neq); %Stiffness matrix
F = zeros(neq, 1); %Force vector
displacements = zeros(size(IDArray)); %u1 and u2 displacements in 1st and 2nd row resp.


%iterate through elements that border the inner radius
for e_1 = 1:n_el1
    e_2 = 1; 
    
    eboundary = n_el1 * (e_2-1) + e_1; %element numbers for inner radius
    
    C_e2 = C_operators(:, :,eboundary);
    
    span1bound = span1vec(e_1,:);
    span2bound = span2vec(e_2, :);
    
    %vector of global node numbers corresponding to local node numbers in
    %element number eboundary
    Avec = IEN(:, eboundary); 
    
    %element NURBS weights and control points
    w_e2 = w(Avec);
    P_ematBoundary = Points(Avec , :);
    
    %element force vector
    Fe = ...
        ElementForceVector2(quad_rule, P_ematBoundary, C_e2, w_e2, span1bound, span2bound); 
    for a = 1: length(Avec)
       A = Avec(a);
       for k = 1:2
          P = ID(k, A); 
          if P == 0
                continue
          end
          p_loc = 2*(a-1) + k;
          F(P) = F(P) + Fe(p_loc); 
       end
    end
    
    
end

 %iterate through elements 
for e = 1:n_el
    %element number of the 1D bezier elements that make up element e
    [e1,e2] = ElementNumber1D(e, n_el1);
    
    %element extraction operator specific to element e
    C_e = C_operators(:,:, e);
    
    span1 = span1vec(e1,:);
    span2 = span2vec(e2, :);
    
    %vector of global node numbers corresponding to local node numbers in
    %element number e
    A_vec = IEN(:, e); 
    
    %element NURBS weights and control points
    w_e = w(A_vec);
    P_e_mat = Points(A_vec , :);
    
    %Element Stiffness matrix
    Ke = ElementStiffness2(quad_rule, P_e_mat, C_e, w_e, span1, span2);
    
    %Assembly
    
    for a = 1: length(A_vec)
        A = A_vec(a);
        for i = 1:2
            P = ID(i, A);
            if P == 0
                continue
            end
            p_loc = 2*(a-1) + i;
            
            for b = 1:length(A_vec)
                B = IEN(b, e);
                for j = 1:2
                Q = ID(j, B);
                if Q == 0
                   continue 
                end
                q_loc = 2*(b-1) + j;
                
                K(P, Q) = K(P, Q) + Ke(p_loc, q_loc);
                end
                
                
            end
        end
        
    end
   
end


%Solve for displacements
d = K\F;


for P = 1:length(d)
    [doflocation, Nodelocation] = find(IDArray == P);
    
    displacements(doflocation, Nodelocation) = d(P);
end

end
