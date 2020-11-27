function node_type_mat = find_node_type(node_type_mat, info_check_vec, in, Nt)
    n = log2(length(info_check_vec));
    nt = log2(Nt);
    
    check_vec = in:in+Nt-1;
    
    % Rate 0 node
    if (all(info_check_vec(check_vec) < 1))
        node_type_mat(n-nt+1,in) = 10;
        
    % Rate 1 node
    elseif(all(info_check_vec(check_vec) > 0))
        node_type_mat(n-nt+1,in) = 20;
    
    % SPC node
    elseif(info_check_vec(check_vec(1)) < 1 && all(info_check_vec(check_vec(2:end)) > 0))
        node_type_mat(n-nt+1,in) = 30;
    
    % Rep node
    elseif(info_check_vec(check_vec(end)) > 0 && all(info_check_vec(check_vec(1:end-1)) < 1))
        node_type_mat(n-nt+1,in) = 40;
    end
    
    if (Nt > 1)
        node_type_mat = find_node_type(node_type_mat, info_check_vec, in, Nt/2);
        node_type_mat = find_node_type(node_type_mat, info_check_vec, in+Nt/2, Nt/2);
    end
end