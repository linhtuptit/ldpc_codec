function p = ldpc_encode_alg(c, H, Z_c)

    c = c(:);
    n = size(H, 2);
    r = size(H, 1);
    k = n - r;
    g = 4 * Z_c;

    % Define sub-matrix for calculation
    A_mat = H(g+1:end, 1:k); 
    B_mat = H(g+1:end, k+1:k+g);
    C_mat = H(1:g, 1:k);
    D_mat = H(1:g, k+1:k+g);

    % Calculate p = [p1, p2]
    D_inv = mod(round(inv(D_mat)), 2);
    Dinv_C = mod(D_inv * C_mat, 2);

    p1 = mod(Dinv_C * c, 2);
    p2 = mod(A_mat * c + B_mat * p1, 2);

    p = [p1; p2];

end