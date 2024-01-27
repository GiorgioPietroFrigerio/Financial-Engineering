% Optimization of the initial conditions for fmincon

% set for the first parameter
first = linspace(0.0005,0.002,200);

% set for the second parameter
second = linspace(0.008,0.009,200);

% options for fmincon
options = optimoptions(@fmincon,'Algorithm','interior-point');

% initialization of matrix of distances
dist_mat = zeros(length(first), length(second));

% construction of the matrix of distances
for i = 1:length(first)
    for j = 1:length(second)
        init_cond = [first(i), second(j)];
        [a, sigma] = fmincon(distance,init_cond,[],[],[],[],[],[],@(parameters) constraints(parameters), options);
        dist = distance([a,sigma]);
        dist_mat(i,j) = dist;
    end
end


% minimization of the matrix of distances
min_dist = min(dist_mat, [], 'all');

% optimal initial conditions
[i, j] = find(dist_mat==min_dist);
init_cond = [first(i), second(j)];

