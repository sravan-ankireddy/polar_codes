function [msg_hat,PM] = decode_SCL(LLR, N, l, info_check_vec, data_pos)
    
    % depth of the polar code tree
    n = log2(N);
    
    % beliefs
    L = zeros(l,n+1,N);
    
    % decisions
    ucap = zeros(l,n+1,N)-2;
    
    % index ordering
    ind_ord_mat = zeros(l,n+1,N);
    
    % vector to check status of node -- left or right or parent propagation
    ns = zeros(1,2*N-1); 
    
    % f function
    f_minsum = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b));
    
    % g function 
    g_minsum = @(a,b,c) b+(1-2*c).*a;
    
    % belief initialisation
    for i = 1:l
        L(i,1,:) = LLR;
    end
    
    % Path Metric
    PM = zeros(l,1);
    
    % start at root
    node = 0; depth = 0;
    
    % stopping criteria
    done = 0; counter = 1;
    % traverse till all bits are decoded
    while (~(done == 1 && depth == -1))
               
        % length of current node
        cur_len = 2^(n-depth);
        
        % position of node in node state vector
        npos = (2^depth-1) + node + 1;
        
        node_type_ind = cur_len*node;
        
        % check for leaf node
        if (depth == n && done == 0)
            % if frozen node
            if (info_check_vec(node+1) == 0)

                ucap(:,n+1,node+1) = 0;
                pen_pos = (L(:,n+1,node+1) < 0);
                PM(pen_pos) = PM(pen_pos) + abs(L(pen_pos,n+1,node+1));
                [PM, ind_ord_l] = sort(PM,'ascend');
                ind_ord_mat(:,n+1,node+1) = ind_ord_l;

            elseif ( counter > log2(l) )
                Lin_temp = [L(:,n+1,node+1);L(:,n+1,node+1)];
                PM_temp = [PM; PM];
                codeword_temp = [zeros(l,1); ones(l,1)];
                pen_pos = ( (Lin_temp < 0) ~= codeword_temp );

                PM_temp(pen_pos) = PM_temp(pen_pos) + abs(Lin_temp(pen_pos));

                [PM_temp, ind_ord] = sort(PM_temp,'ascend');

                ucap(:,n+1,node+1) = codeword_temp(ind_ord(1:l));

                PM = PM_temp(1:l); 
                ind_ord = ind_ord(1:l);
                ind_ord(ind_ord > l) = ind_ord(ind_ord > l) - l;
                ind_ord_mat(:,n+1,node+1) = ind_ord;
                counter = counter + 1;
                
            elseif (counter == 1 || counter == 2 || counter == 3 || counter == 4)
                
                if (counter == 1)
                    codeword = [zeros(l/2,1); ones(l/2,1)];
                elseif (counter == 2)
                    codeword = [zeros(l/4,1); ones(l/4,1)];
                    codeword = repmat(codeword,2,1);
                elseif (counter == 3)
                    codeword = [zeros(l/8,1); ones(l/8,1)];
                    codeword = repmat(codeword,4,1);
                elseif (counter == 4)
                    codeword = [zeros(l/16,1); ones(l/16,1)];
                    codeword = repmat(codeword,8,1);
                end
                    
                ucap(:,n+1,node+1) = codeword;
                Lin = squeeze(L(:,n+1,node+1));
                pen_pos = ( (Lin < 0) ~= codeword );
                PM(pen_pos) = PM(pen_pos) + abs(Lin(pen_pos));
                [PM, ind_ord_l] = sort(PM,'ascend');
                ind_ord_mat(:,n+1,node+1) = ind_ord_l;
                ucap(:,n+1,node+1) = codeword(ind_ord_l);
                counter = counter + 1;
            end
            % check for last leaf node
            if node == (N-1)
                done = 1; node = floor(node/2); depth = depth - 1;

                % move back to parent node
            else
                node = floor(node/2); depth = depth - 1;
            end
            
        % non-leaf node
        else
            
            % propagate to left child
            if ns(npos) == 0
                
                % incoming beliefs
                Ln = L(:,depth+1,node_type_ind + 1:node_type_ind + cur_len);
                
                % next node: left child
                node = 2*node; depth = depth + 1; 
                
                % incoming belief length for left child
                cur_len = floor(cur_len / 2);
                
                % calculate and store LLRs for left child
                L(:,depth+1,node_type_ind+1:node_type_ind+cur_len) = f_minsum(Ln(:,1,1:cur_len),Ln(:,1,cur_len+1:end));
                
                % mark as left child visited
                ns(npos) = 1;
            else
                % propagate to right child
                if ns(npos) == 1
                    
                    % incoming beliefs
                    ind_cur_ord = ind_ord_mat(:,depth+2,node_type_ind + 1);
                    Ln = L(ind_cur_ord,depth+1,node_type_ind+1:node_type_ind+cur_len);
                    
                    % left child
                    lnode = 2*node; ldepth = depth + 1; 
                    ltemp = cur_len/2;
                    
                    % incoming decisions from left child
                    ucapn = ucap(:,ldepth+1,ltemp*lnode+1:ltemp*(lnode+1));
                    
                    % next node: right child
                    node = node *2 + 1; depth = depth + 1;
                    
                    % incoming belief length for right child
                    cur_len = floor(cur_len / 2);
                    
                    % calculate and store LLRs for right child
                    L(:,depth+1,cur_len*node+1:cur_len*(node+1)) = g_minsum(Ln(:,1,1:cur_len),Ln(:,1,cur_len+1:end),ucapn);
                    
                    % mark as right child visited
                    ns(npos) = 2;
                    
                % calculate beta propagate to parent node
                else
                    % ordering 
                    ind_ord_temp = ind_ord_mat(:, depth+2, node_type_ind + 1);
                    
                    ind_ord_temp2 = ind_ord_mat(:, depth+2, node_type_ind + 1 + floor(cur_len/2));
                                        
                    ind_ord_mat(:,depth+1,node_type_ind + 1) = ind_ord_temp(ind_ord_temp2);
                    
                    % left and right child
                    lnode = 2*node; rnode = 2*node + 1; cdepth = depth + 1;
                    ctemp = cur_len/2;
                    
                    % incoming decisions from left child
                    i_temp_ord = ind_ord_mat(:,cdepth+1,node_type_ind + 1 + cur_len/2);
                    ucapl = ucap(:,cdepth+1,ctemp*lnode+1:ctemp*(lnode+1));
                    
                    % incoming decisions from right child
                    ucapr = ucap(:,cdepth+1,ctemp*rnode+1:ctemp*(rnode+1));
                    
                    % combine
                    ucap(:,depth+1,node_type_ind+1:node_type_ind+cur_len) = cat(3, mod(ucapl(i_temp_ord,:,:)+ucapr,2), ucapr);
                    
                    % update to index of parent node
                    node = floor(node/2); depth = depth - 1;
                end
            end
        end
    end
    
    if l == 1
        msg_hat = squeeze(ucap(:,n+1,data_pos))';
    else
        msg_hat = squeeze(ucap(:,n+1,data_pos));
    end

end