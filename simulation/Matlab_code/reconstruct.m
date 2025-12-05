function A = reconstruct(T, X, M, W, opts)
%% input
% T: ground-truth tensor (dim: ng-nx-ny-nz)
% X: cell array of gene expression matrices 
%    along spatial coordinates
% M: mask tensor (dim: 1-nx-ny-nz)
% W: cell array of graphs (PPI graph and 
%    chain graphs along spatial coordinates)  
% opts: parameters
%    - rank: tensor CPD rank
%    - alpha: hyperparameter
%    - beta: hyperparameter
%    - lambda: hyperparameter
%    - stopcrit: stop criteria
%    - maxiters: maximum iteration
%% output
% A: cell array of CPD component matrices

eps = 1e-10; % numerical stability

n = size(T); % tensor dimensions
n_mode = length(n); % tensor modes

% normalize graph
D = cell(n_mode, 1);
for i = 1:n_mode
    W{i} = W{i} - diag(diag(W{i}));
    d = sum(W{i}, 2);
    d(d ~= 0) = (d(d ~= 0)) .^ -(0.5);
    W{i} = W{i} .* d;
    W{i} = d' .* W{i};
    D{i} = diag(sum(W{i}, 2));
end

% initialize tensor CPD components
% and auxiliary variables
A = cell(1, n_mode);
ATA = cell(1, n_mode);
aaT = cell(1, n_mode);
WA = cell(1, n_mode);
DA = cell(1, n_mode);
ATWA = cell(1, n_mode);
ATDA = cell(1, n_mode);
for i = 1:n_mode
    rng(0);
    A{i} = rand(n(i), opts.rank);
    ATA{i} = A{i}' * A{i};
    aaT{i} = sum(A{i}, 1)' * sum(A{i}, 1);
    WA{i} = W{i} * A{i};
    DA{i} = diag(D{i}) .* A{i};
    ATWA{i} = A{i}' * WA{i};
    ATDA{i} = A{i}' * DA{i};
end

% convert expression matrices to tensors
X_tensor = cell(1, n_mode);
for i = 2:n_mode
    [rows, cols, vals] = find(X{i});
    subs = ones(length(rows), n_mode);
    subs(:, [1, i]) = [rows, cols];
    dims = ones(1, n_mode);
    dims(1) = n(1);
    dims(i) = n(i);
    X_tensor{i} = sptensor(subs, vals, dims);
end

% % collapse tensor components
% A_collap_all = cell(1, n_mode);
% for i = 1:n_mode
%     A_collap_all{i} = sum(A{i}, 1);
% end

% compute scaling factors
mode_idx = 1:n_mode;
fac = cell(n_mode, n_mode);
for i = 1:n_mode-1
    for j = i+1:n_mode
        fac{i, j} = prod(n(setdiff(mode_idx, [i, j])));
    end
end


% core model
for iter = 1:opts.maxiters

    Aold = A;

    % update gene component
    num = zeros(size(A{1}));
    denom = zeros(size(A{1}));
    % % NEW ADDED
    A_collap_all = collapse_all(A, n_mode);
    % % 

    % update J1 negative
    theta = cell(1, n_mode);
    J1n = zeros(size(A{1}));
    
    for i = 2:n_mode
        A_collap = A_collap_all;
        A_collap{1} = A{1};
        A_collap{i} = A{i};
        theta{i} = mttkrp(X_tensor{i}, A_collap, 1);
        J1n = J1n + theta{i}/fac{1, i};
    end
    num = num + J1n + eps;

    % update J1 positive
    phi = cell(1, n_mode);
    phi_sum = zeros(opts.rank, opts.rank);
    J1p = zeros(size(A{1}));
    for i = 2:n_mode
        phi{i} = ATA{i};
        for j = 2:n_mode
            if j ~= i
                phi{i} = phi{i} .* aaT{j};
            end
        end
        phi_sum = phi_sum + phi{i}/(fac{1, i}^2);
    end
    J1p = J1p + A{1} * phi_sum;
    denom = denom + J1p + eps;

    % update J2 positive
    A_collap = A;
    % % A_collap{1} = sum(A_collap{1}, 1);
    A_collap{1} = A_collap_all{1};
    J2p = mttkrp(ktensor(A_collap) .* M, A_collap, 1);
    J2p = repmat(J2p, n(1), 1);
    denom = denom + double(opts.alpha) * J2p/double((n(1)^2)) + eps;

    % update J3 negative
    psi = cell(1, n_mode);
    psi_sum = zeros(opts.rank, opts.rank);
    ATA_prod = ones(opts.rank, opts.rank);
    for i = 2:n_mode
        ATA_prod = ATA_prod .* ATA{i};
        psi{i} = ATWA{i};
        for j = 2:n_mode
            if j~= i
                psi{i} = psi{i} .* ATA{j};
            end
        end
        psi_sum = psi_sum + psi{i};
    end
    J3n =  WA{1} * ATA_prod + A{1} * psi_sum;
    num = num + double(opts.beta) * J3n + eps;

    % update J3 positive
    psi = cell(1, n_mode);
    psi_sum = zeros(opts.rank, opts.rank);
    ATA_prod = ones(opts.rank, opts.rank);
    for i = 2:n_mode
        ATA_prod = ATA_prod .* ATA{i};
        psi{i} = ATDA{i};
        for j = 2:n_mode
            if j~= i
                psi{i} = psi{i} .* ATA{j};
            end
        end
        psi_sum = psi_sum + psi{i};
    end
    J3p =  DA{1} * ATA_prod + A{1} * psi_sum;
    denom = denom + double(opts.beta) * J3p + eps;

    % update J4 positive
    J4p = A{1};
    denom = denom + double(opts.lambda) * J4p + eps;
    
    % calculate new gene component
    A{1} = A{1} .* (num ./ denom);
    
    % update auxiliary variables
    ATA{1} = A{1}' * A{1};
    aaT{1} = sum(A{1}, 1)' * sum(A{1}, 1);
    WA{1} = W{1} * A{1};
    DA{1} = diag(D{1}) .* A{1};
    ATWA{1} = A{1}' * WA{1};
    ATDA{1} = A{1}' * DA{1};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update spatial component

    for i = 2:n_mode

        num = zeros(size(A{i}));
        denom = zeros(size(A{i}));
        % % NEW ADDED 
        A_collap_all = collapse_all(A, n_mode);
        % %

        % update J1 negative
        theta = cell(1, n_mode);
        J1n = zeros(size(A{i}));
        A_collap = A_collap_all;
        A_collap{1} = A{1};
        A_collap{i} = A{i};
        theta{1} = mttkrp(X_tensor{i}, A_collap, i);
        J1n = J1n + theta{1}/fac{1, i};

        for j = 2:n_mode
            if j ~= i
                A_collap = A_collap_all;
                A_collap{1} = A{1};
                A_collap{j} = A{j};
                theta{j} = mttkrp(X_tensor{j}, A_collap, i);
                theta{j} = repmat(theta{j}, n(i), 1);
                J1n = J1n + theta{j}/fac{1, j};
            end
        end
        num = num + J1n + eps;

        % update J1 positive
        phi = cell(1, n_mode);
        J1p = zeros(size(A{i}));
        phi{1} = ATA{1};
        for j = 2:n_mode
            if j~=i
                phi{1} = phi{1} .* aaT{j};
            end
        end
        J1p = J1p + A{i} * (phi{1}/(fac{1, i}^2));

        phi_sum = zeros(opts.rank, opts.rank);
        for j = 2:n_mode
            if j ~= i
                phi{j} = ATA{1} .*  ATA{j};
                for k = 2:n_mode
                    if k~=i && k~=j
                        phi{j} = phi{j} .* aaT{k};
                    end
                end
                phi_sum = phi_sum + phi{j}/(fac{1, j}^2);
            end
        end
        J1p = J1p + repmat(sum(A{i}, 1) * phi_sum, n(i), 1);
        denom = denom + J1p + eps;

        % update J2 positive
        A_collap = A;
        A_collap{1} = A_collap_all{1};
        J2p = mttkrp(ktensor(A_collap) .* M, A_collap, i);
        denom = denom + double(opts.alpha) * J2p/double((n(1)^2)) + eps;

        % update J3 negative
        psi = cell(1, n_mode);
        psi_sum = zeros(opts.rank, opts.rank);
        ATA_prod = ones(opts.rank, opts.rank);
        for j = 1:n_mode
            if j~=i
                ATA_prod = ATA_prod .* ATA{j};
                psi{j} = ATWA{j};
                for k = 1:n_mode
                    if k ~= i && k ~= j
                        psi{j} = psi{j} .* ATA{k};
                    end
                end
                psi_sum = psi_sum + psi{j};
            end
        end
        J3n =  WA{i} * ATA_prod + A{i} * psi_sum;
        num = num + double(opts.beta) * J3n + eps;

        % update J3 positive
        psi = cell(1, n_mode);
        psi_sum = zeros(opts.rank, opts.rank);
        ATA_prod = ones(opts.rank, opts.rank);
        for j = 1:n_mode
            if j~=i
                ATA_prod = ATA_prod .* ATA{j};
                psi{j} = ATDA{j};
                for k = 1:n_mode
                    if k ~= i && k ~= j
                        psi{j} = psi{j} .* ATA{k};
                    end
                end
                psi_sum = psi_sum + psi{j};
            end
        end
        J3p =  DA{i} * ATA_prod + A{i} * psi_sum;
        denom = denom + double(opts.beta) * J3p + eps;

        % update J4 positive
        J4p = A{i};
        denom = denom + double(opts.lambda) * J4p + eps;
        
        % calculate new spatial component
        A{i} = A{i} .* (num ./ denom);
        
        % update auxiliary variables
        ATA{i} = A{i}' * A{i};
        aaT{i} = sum(A{i}, 1)' * sum(A{i}, 1);
        WA{i} = W{i} * A{i};
        DA{i} = diag(D{i}) .* A{i};
        ATWA{i} = A{i}' * WA{i};
        ATDA{i} = A{i}' * DA{i};

    end

    res = compute_res(A, Aold);
    disp(['training...residual: ', num2str(res)]);
    disp(['training...iteration: ', num2str(iter)]);
    if iter >= opts.miniters && res < opts.stopcrit
        break;
    end

end

end

% Auxiliary function
function A_collap_all = collapse_all(A, n_mode)
% collapse tensor components
A_collap_all = cell(1, n_mode);
for i = 1:n_mode
    A_collap_all{i} = sum(A{i}, 1);
end
end

function res = compute_res(A, Aold)
% compute residual
eps = 1e-10; % numerical stability
res_num = 0;
res_denom = eps;
for i = 1:length(A)
    res_num = res_num + sum(sum((A{i} - Aold{i}).^2));
    res_denom = res_denom + sum(sum(Aold{i}.^2));
end
res = sqrt(res_num/res_denom);
end
