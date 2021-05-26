function wd = ldpc_decode_alg(w, H, max_iter)

    V_temp = zeros(1, length(w));
    CV_mesg = zeros(size(H, 1), size(H, 2));
    VC_mesg = zeros(size(H, 2), size(H, 1));
    V_out = zeros(1, length(w));

    % Initialization Processing ...
    for idx = 1:1:length(w)
        V_temp(idx) = w(idx);
    end
    clear idx;

    for idx = 1:1:size(VC_mesg, 1)
        for jdx = 1:1:size(VC_mesg, 2)
            if H(jdx, idx) == 1
                VC_mesg(idx, jdx) = V_temp(idx);
            end
        end
    end
    clear idx; clear jdx;

    % START FOR LOOP ...
    for iter = 1:1:max_iter
    % Check to Variable node processing ...
        for idx = 1:1:size(CV_mesg, 1)
            for jdx = 1:1:size(CV_mesg, 2)
                if H(idx, jdx)
                    min_cv = 10000;
                    prod_cv = 1;
                    for kdx = 1:1:size(CV_mesg, 2)
                        if H(idx, kdx) && kdx ~= jdx
                            if abs(VC_mesg(kdx, idx)) < min_cv
                                min_cv = abs(VC_mesg(kdx, idx));
                            end
                            prod_cv = prod_cv*sign(VC_mesg(kdx, idx));
                        end
                    end
                    CV_mesg(idx, jdx) = prod_cv * min_cv;
                end
            end
        end
        clear idx; clear jdx; clear kdx;
        % 
        % Variable to Check node processing ...
        for idx = 1:1:size(VC_mesg, 1)
            for jdx = 1:1:size(VC_mesg, 2)
                if H(jdx, idx)
                    sum_vc = 0;
                    for kdx = 1:1:size(VC_mesg, 2)
                        if H(kdx, idx) && jdx ~= kdx
                            sum_vc = sum_vc + CV_mesg(kdx,idx);
                        end
                    end
                    VC_mesg(idx, jdx) = V_temp(idx) + sum_vc;
                end
            end
        end

        % Variable node updating ...
        for idx = 1:1:length(w)
            V_out(idx) = V_temp(idx) + sum(CV_mesg(:,idx));
        end
        clear idx;

        % Hard Decode and Syndrome Checking ...
        wd = hard_func(V_out);
        CkSy = mod(H*wd', 2);
        if sum(CkSy) == 0 || iter == max_iter
            fprintf("\n Number of iterations are: %d \n", iter);
            break;
        end
        clear idx;
    % END FOR LOOP
    end
    
    function out = hard_func(in)
        out = zeros(1, length(in));
        for tdx = 1:1:length(in)
            if in(tdx) >= 0; out(tdx) = 0;
            elseif in(tdx) < 0; out(tdx) = 1;
            end
        end
    end

end