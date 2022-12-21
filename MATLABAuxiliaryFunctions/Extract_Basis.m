function [n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2)
% SPECIAL THANKS TO MINYEONG KIM FOR ALLOWING ME TO RE-USE HIS CODE

% Function to construct the element extraction operator Ce and IEN array 
% for 2D where, 
% p_1, p_2 : the polynomial degree of each dir.
% n_1, n_2 : the number of basis in each dir.
% Xi_1, Xi_2 : the univariate knot vector
% n_el : the number of Bezier elements
% C_operators : array storing the ELEMENT extraction operator
% IEN : array mapping local basis function numbers to
% global basis function number

% 1. PROCEED BEZIER EXTRACTION IN EACH DIRECTION
[D1.n_el,D1.C_operators,D1.IEN] = Extract_Basis_1D(p_1,n_1,Xi_1);
[D2.n_el,D2.C_operators,D2.IEN] = Extract_Basis_1D(p_2,n_2,Xi_2);

% 3. DETERMINE TOTAL BASIS NUMBER n_el FOR GLOBAL SPACE
n_el = D1.n_el*D2.n_el;

% 2. PERFORM KRONECKER TENSOR PRODUCT BETWEEN EXTRACTION OPERATORS
C_operators = zeros((p_2+1)*(p_1+1), (p_2+1)*(p_1+1), n_el);
for e2=1:D2.n_el
    for e1=1:D1.n_el
        C_operators(:,:,D1.n_el*(e2-1)+e1) = ...
        kron(D2.C_operators(:,:,e2), D1.C_operators(:,:,e1));
    end
end



% 4. DETERMINE IEN ARRAY MAPS
% It defines the maps from local basis # and element # onto global
% basis #. Column and row of the array denotes element # and local
% basis # respectively. Here, indexing: A = n_1*(j-1)+i
IENj = kron(n_1*(D2.IEN-1), zeros(p_1+1, D1.n_el)+1);
IENi = kron(zeros(p_2+1, D2.n_el)+1, D1.IEN);
IEN = IENi+IENj;


%--------------------------------------------------------------------------
    % Supplementary Function Body (1D Bezier Extraction)
    function [n_el,C_operators,IEN] = Extract_Basis_1D(p,n,Xi)
        % 1. INITIALIZE FOR BEZIER DECOMPOSITION
        % Here, knots will be inserted so that
        % each internal knot has multiplicity of p.
        Xii = Xi; % Store input Xi temporarily
        Xilist = unique(Xii); %
        j=1;
        i=p+2;
        mult(1,j) = 1;
        % Compute multiplicity for each internal knot and stored in "mult"
        % Also, specifies the index of knot vector where internal knot changes
        % value in "loc"
        while i>=p+2 && i<length(Xi)-p
            if Xi(i)==Xi(i+1)
                mult(1,j) = mult(1,j)+1;
                i = i+1;
            elseif Xi(i)<Xi(i+1)
                loc(1,j) = i;
                i = i+1;
                j = j+1;
                mult(1,j) = 1;
            end   
        end
        mult(j) = [];
        
        amult = p*ones(1,length(mult))-mult; % necessary # of multiplicity
        % 2. DETERMINE KNOT VECTOR TO ADDED
        add_Xi = []; % initialize add_Xi
        for i=1:length(amult)
            add_Xi = [add_Xi Xii(loc(i))*ones(1,amult(i))];
        end
        
        % 3. KNOT INSERTION ALGORITHM:
        % Retrived from the function, NURBS_Curve_Refine
        % Initialize the new knot vector array
        new_Xi = Xi;
        %len_Xi = length(Xi);
        C = 1;
        for i=1:length(add_Xi)
            for j=1:length(new_Xi)
                if add_Xi(i)>=new_Xi(j) && add_Xi(i)<new_Xi(j+1)
                % Knot insertion and save to the temporal directory
                b = [new_Xi(1:j) add_Xi(i) new_Xi(j+1:length(new_Xi))];
                % and store its location.
                k = j;
                end
            end
            % If the process ends with i-th element of add_Xi, store array in
            % temporal directory to the new knot vector array.
            old_Xi = new_Xi; % for determining the alpha(j) below;
            new_Xi = b;
            % Determine the no. of basis ftn’s
            new_n = length(new_Xi)-p-1;
        
            % Determine the set of alpha in this case
            alpha = zeros(1, new_n);
            for j=1:new_n
                if j>=1 && j<k-p+1
                    alpha(j) = 1;
                elseif j>=k-p+1 && j<k+1
                    alpha(j) = (add_Xi(i)-old_Xi(j))/(old_Xi(j+p)-old_Xi(j));
                elseif j>=k+1 && j<=new_n
                    alpha(j) = 0;
                end
            end
            
             % 4. CONSTRUCT BEZIER EXTRACTION OPERATOR AT j-TH KNOT INSERTION
             Cj1 = [alpha(1:new_n-1).*eye(new_n-1) zeros(new_n-1,1)];
             Cj2 = [zeros(new_n-1,1) eye(new_n-1)-alpha(2:new_n).*eye(new_n-1)];
             Cj = Cj1+Cj2; % (n+j-1)*(n+j) matrix
 
            C = C*Cj; % implicitly calculate the extraction operator
            clear Cj1 Cj2 Cj   
        end
        
      
        if isempty(add_Xi) == 1
            % If there is no knot insertion, let C be:
            C = eye(n, n);
        end

        % 5. DETERMINE THE NUMBER OF ELEMENT
        n_el = length(mult)+1;
        
        % 6. DETERMINE THE KNOTS SUPPORTING EACH B-SPLINE CURVE
        % Additional code for constructing element extract operator for
        % non-uniform knot vector.
        supp = zeros(n, p+2);
        for i=1:n
            supp(i,:) = Xii(i:i+p+1);
        end
        % 7. DECOMPOSING THE BEZIER EXTRACTION OPERATOR INTO C_operators(:,:,e)
        C_operators = zeros(p+1, p+1 ,n_el);
        IEN = zeros(p+1, n_el);
        for e=1:n_el
            j = 1;
            % Determine the overlaping basis for each Bezier element
            while j<=n && sum(supp(j,:)==Xilist(e+1))==0
                j =j+1;
            end
            % Determine the element extraction operators
            C_operators(:,:,e) = C(j:j+p,p*(e-1)+1:p*e+1);
            
            
            % 8. DETERMINE THE IEN(a,e) ARRAY
            for i=1:p+1
            IEN(i,e) = j+i-1;
            end
            %aloc(e,:) = [j e]
        end

    end

end

