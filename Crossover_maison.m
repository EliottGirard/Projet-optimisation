function xoverKids = Crossover_maison(parents,options,GenomeLength,~,~,thisPopulation, n1)
% Fonction crossover fait maison, garde C et i ensembles, laplace crossover
% Call : xoverKids = myfun(parents, options, nvars, FitnessFcn, ...,unused,thisPopulation)

% Modification de crossoverscattered

% Données pour crossover laplacien
a_lapl = 1;
b_lapl = 0.5;

% How many children to produce?
nKids = length(parents)/2;
% Extract information about linear constraints, if any
linCon = options.LinearConstr;
constr = ~isequal(linCon.type,'unconstrained');
% Allocate space for the kids
xoverKids = zeros(length(parents),GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;

% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Randomly select half of the genes from each parent
    % ------> Modif : Create 2 kids based on each parents
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    for j = 1:n1
        i1 = double(thisPopulation(r1,j+n1) > 0.5);
        i2 = double(thisPopulation(r2,j+n1) > 0.5);
        
        xoverKids(i,j+n1) = i1;
        xoverKids(i+nKids,j+n1) = i2;
        
        if(i1&&i2)
            beta_lapl = a_lapl + (2*(rand>0.5) - 1)*b_lapl*(rand>0.5);
            xoverKids(i, j) = thisPopulation(r1, j) + beta_lapl*abs(thisPopulation(r1, j) - thisPopulation(r2, j));
            xoverKids(i+nKids, j) = thisPopulation(r2, j) + beta_lapl*abs(thisPopulation(r1, j) - thisPopulation(r2, j));
        else
            xoverKids(i, j) = i1*thisPopulation(r1, j);
            xoverKids(i+nKids, j) = i2*thisPopulation(r2, j);
        end
    end
    % Pas de test de faisabilité
end
end

