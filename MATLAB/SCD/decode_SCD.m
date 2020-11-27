function msg_hat = decode_SCD(LLR, N, info_check_vec, data_pos)
    
    % depth of the polar code tree
    n = log2(N);
    
    % beliefs
    L = zeros(n+1,N);
    
    % decisions
    ucap = zeros(n+1,N);
    
    % vector to check status of node -- left or right or parent propagation
    ns = zeros(1,2*N-1); 
    
    % f function
    f_minsum = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b));
    
    % g function 
    g_minsum = @(a,b,c) b+(1-2*c).*a;
    
    % belief initialisation
    L(1,:) = LLR;
    
    % start at root
    node = 0; depth = 0;
    
    % stopping criteria
    done = 0;
    % traverse till all bits are decoded
    while (done == 0) 
        
        % check for leaf node
        if depth == n
            % if frozen node
            if (info_check_vec(node+1) == 0)
                ucap(n+1,node+1) = 0;
            else
                if (L(n+1,node+1) >= 0)
                    ucap(n+1,node+1) = 0;
                else
                    ucap(n+1,node+1) = 1;
                end
            end
            
            % check for last leaf node
            if node == (N-1)
                done = 1;
            % move back to parent node
            else
                node = floor(node/2); depth = depth - 1;
            end
            
        % non-leaf node
        else
            % position of node in node state vector
            npos = (2^depth-1) + node + 1;
            
            % propagate to left child
            if ns(npos) == 0
                % length of current node
                cur_len = 2^(n-depth);
                % incoming beliefs
                Ln = L(depth+1,cur_len*node+1:cur_len*(node+1));
                % next node: left child
                node = 2*node; depth = depth + 1; 
                % incoming belief length for left child
                cur_len = cur_len / 2;
                % calculate and store LLRs for left child
                L(depth+1,cur_len*node+1:cur_len*(node+1)) = f_minsum(Ln(1:cur_len),Ln(cur_len+1:end));
                % mark as left child visited
                ns(npos) = 1;
            else
                % propagate to right child
                if ns(npos) == 1
                    % length of current node
                    cur_len = 2^(n-depth);
                    % incoming beliefs
                    Ln = L(depth+1,cur_len*node+1:cur_len*(node+1));
                    % left child
                    lnode = 2*node; ldepth = depth + 1; 
                    ltemp = cur_len/2;
                    % incoming decisions from left child
                    ucapn = ucap(ldepth+1,ltemp*lnode+1:ltemp*(lnode+1));
                    % next node: right child
                    node = node *2 + 1; depth = depth + 1;
                    % incoming belief length for right child
                    cur_len = cur_len / 2;
                    % calculate and store LLRs for right child
                    L(depth+1,cur_len*node+1:cur_len*(node+1)) = g_minsum(Ln(1:cur_len),Ln(cur_len+1:end),ucapn);
                    % mark as right child visited
                    ns(npos) = 2;
                % calculate beta propagate to parent node
                else
                    % length of current node
                    cur_len = 2^(n-depth);
                    % left and right child
                    lnode = 2*node; rnode = 2*node + 1; cdepth = depth + 1;
                    ctemp = cur_len/2;
                    % incoming decisions from left child
                    ucapl = ucap(cdepth+1,ctemp*lnode+1:ctemp*(lnode+1));
                    % incoming decisions from right child
                    ucapr = ucap(cdepth+1,ctemp*rnode+1:ctemp*(rnode+1));
                    % combine
                    ucap(depth+1,cur_len*node+1:cur_len*(node+1)) = [mod(ucapl+ucapr,2) ucapr];
                    % update to index of parent node
                    node = floor(node/2); depth = depth - 1;
                end
            end
        end
    end
    
    msg_hat = ucap(n+1,data_pos);
    
end