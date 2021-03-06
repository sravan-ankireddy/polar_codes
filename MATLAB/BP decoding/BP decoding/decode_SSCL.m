function [y_hat, ind_ord, PM, counter, beta] = decode_SSCL(LLR, PM, in, counter, y_hat, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array)
    % list size
    l = size(LLR,1);
    
    % Length of current sub-vector
	N = size(LLR,2);
    
    % currect subvec ind
    ind_vec = in:in+N-1;
    
    ind_ord = 1:2*l;
    if ( N == 1 )
                
        if ( Rate_0_node_pos_array(in, N) > 0 )
            codeword = zeros(l,1);
            for i_list = 1:l
                if (codeword(i_list) ~= 0.5*(1 - sign(LLR(i_list))))
                    PM(i_list) = PM(i_list) + abs(LLR(i_list));
                end
            end
            
            y_hat(:,in) = codeword;
            [PM, ind_PM] = sort(PM, 'ascend');
            ind_ord = ind_PM;
            ind_ord = [ind_ord; ind_ord];
            y_hat = y_hat(ind_ord(1:l),:);
            beta = y_hat(:,in);
            
        else
            
            if (counter > log2(l))
                
                codeword_temp = zeros(2*l,1);
                y_hat_temp = [y_hat; y_hat];
                LLR_temp = [LLR; LLR];
                PM_temp = [PM; PM];

                for i_list = 1:l
                    codeword_temp(i_list) = 0;
                    if (codeword_temp(i_list) ~= 0.5*(1 - sign(LLR_temp(i_list))))
                        PM_temp(i_list) = PM_temp(i_list) + abs(LLR_temp(i_list));
                    end
                end
                
                % Path 1
                for i_list = l+1:2*l
                    codeword_temp(i_list) = 1;
                    if (codeword_temp(i_list) ~= 0.5*(1 - sign(LLR_temp(i_list))))
                        PM_temp(i_list) = PM_temp(i_list) + abs(LLR_temp(i_list));
                    end
                end
                
                % Pruning
                [PM, ind_PM] = sort(PM_temp, 'ascend');
                PM = PM(1:l);
                ind_ord = ind_PM;
                
                y_hat_temp(:,in) = codeword_temp;
                y_hat = y_hat_temp(ind_PM(1:l),:);
                counter = counter + 1;
                beta = y_hat(:,in);
                        
            elseif (counter == 1)
                                
                codeword = zeros(l,1);
                for i_codeword = l/2+1:l
                    codeword(i_codeword) = 1;
                end
                
                for i_list = 1:l
                    if (codeword(i_list) ~= 0.5*(1 - sign(LLR(i_list))))
                        PM(i_list) = PM(i_list) + abs(LLR(i_list));
                    end
                end
                
                counter = counter + 1;
                y_hat(:,in) = codeword;
                [PM, ind_PM] = sort(PM, 'ascend');
                ind_ord = ind_PM;
                ind_ord = [ind_ord; ind_ord];
                y_hat = y_hat(ind_ord(1:l),:);
                beta = y_hat(:,in);
                
            elseif (counter == 2)
                                
                codeword = zeros(l,1);
                for i_codeword = 1:l
                    if (mod(i_codeword,2) == 0)
                        codeword(i_codeword) = 1;
                    end
                end
                
                for i_list = 1:l
                    if (codeword(i_list) ~= 0.5*(1 - sign(LLR(i_list))))
                        PM(i_list) = PM(i_list) + abs(LLR(i_list));
                    end
                end
                
                counter = counter + 1;
                y_hat(:,in) = codeword;
                [PM, ind_PM] = sort(PM, 'ascend');
                ind_ord = ind_PM;
                ind_ord = [ind_ord; ind_ord];
                y_hat = y_hat(ind_ord(1:l),:);
                beta = y_hat(:,in);
                
            end
        end
        
    elseif ( Rate_0_node_pos_array(in, length(ind_vec)) > 0 )
        codeword = zeros(l,N);
        for i_list = 1:l
            PM(i_list) = PM(i_list) + sum(log(1 + exp(-LLR(i_list,:))));
%             PM(i_list) = PM(i_list) + 0.5 * sum( sign(LLR(i_list,:).*LLR(i_list,:) - LLR(i_list,:)) );
        end
        
        y_hat(:,in:in+N-1) = codeword;
        [PM, ind_PM] = sort(PM, 'ascend');
        ind_ord = ind_PM;
        ind_ord = [ind_ord; ind_ord];
        y_hat = y_hat(ind_ord(1:l),:);
        beta = y_hat(:,in:in+N-1);
        
    elseif ( Rep_node_pos_array(in, length(ind_vec)) > 0 )
        if (counter > log2(l))
            codeword_temp = zeros(2*l,N);
            codeword_temp(l+1:2*l,:) = ones(l,N);
            LLR_temp = [LLR; LLR];
            PM_temp = [PM; PM];

            PM_temp = PM_temp + sum(log(1 + exp(-(1 - 2*codeword_temp).*LLR_temp)), 2);

            % Pruning
            [PM, ind_PM] = sort(PM_temp, 'ascend');
            PM = PM(1:l);
            ind_ord = ind_PM;
            
            y_hat_temp = [y_hat; y_hat];
            y_hat_temp(:,in:in+N-1) = codeword_temp;
            y_hat = y_hat_temp(ind_PM(1:l),:);
            counter = counter + 1;
            beta = y_hat(:,in:in+N-1);
            for i_en = 1:l
                y_hat(i_en,in:in+N-1) = encode(y_hat(i_en,in:in+N-1));
            end
        elseif (counter == 1)
            codeword = zeros(l,N);
            codeword(l/2+1:l,:) = ones(l/2,N);
            for i_list = 1:l
                PM(i_list) = PM(i_list) + sum(log(1 + exp(-(1 - 2*codeword(i_list,:)).*LLR(i_list,:))));
            end

            [PM, ind_PM] = sort(PM, 'ascend');
            ind_ord = [ind_PM; ind_PM];
            counter = counter + 1;
            y_hat(:,in:in+N-1)= codeword;
            y_hat = y_hat(ind_ord(1:l),:);
            beta = y_hat(:,in:in+N-1);
            for i_en = 1:l
                y_hat(i_en,in:in+N-1) = encode(y_hat(i_en,in:in+N-1));
            end
        elseif (counter == 2)
            codeword = zeros(l,N);
            for i_codeword = 1:l
                if (mod(i_codeword,2) == 0)
                    codeword(i_codeword,:) = ones(1,N);
                end
            end
            for i_list = 1:l
                PM(i_list) = PM(i_list) + sum(log(1 + exp(-(1 - 2*codeword(i_list,:)).*LLR(i_list,:))));
            end
            
            [PM, ind_PM] = sort(PM, 'ascend');
            ind_ord = [ind_PM; ind_PM];
            counter = counter + 1;
            y_hat(:,in:in+N-1)= codeword;
            y_hat = y_hat(ind_ord(1:l),:);
            beta = y_hat(:,in:in+N-1);
            for i_en = 1:l
                y_hat(i_en,in:in+N-1) = encode(y_hat(i_en,in:in+N-1));
            end
        end
        
    elseif ( Rate_1_node_pos_array(in, length(ind_vec)) > 0 )
        codeword_temp = zeros(2*l,N) - 2;
        LLR_temp = [LLR; LLR];
        PM_temp = [PM; PM];
        y_hat_temp = [y_hat; y_hat];
        
        codeword = zeros(2*l,1);
        codeword(l+1:2*l,:) = ones(l,1);
        ind_ord = 1:2*l;
        for i_N = 1:N

            codeword_temp(:,i_N) = codeword;
            PM_temp = PM_temp + log(1 + exp( -(1 - 2*codeword_temp(:,i_N)).*LLR_temp(:,i_N)));
            
             % Pruning
            [PM, ind_PM] = sort(PM_temp, 'ascend');
            PM = PM(1:l);
            y_hat_temp(:,in+i_N-1) = codeword;
            y_hat = y_hat_temp(ind_PM(1:l),:);
            LLR = LLR_temp(ind_PM(1:l),:);
            codeword_temp = [codeword_temp(ind_PM(1:l),:); codeword_temp(ind_PM(1:l),:)];
            ind_ord = ind_ord([ind_PM(1:l) ind_PM(1:l)]);

            counter = counter + 1;
            LLR_temp = [LLR; LLR];
            PM_temp = [PM; PM];
            y_hat_temp = [y_hat; y_hat];

        end
        beta = y_hat(:,in:in+N-1);
        for i_en = 1:l
            y_hat(i_en,in:in+N-1) = encode(beta(i_en,:));
        end
        
    else
        
        % Left fork
        LLR_left = zeros(l,N/2);
        for i_left = 0:N/2 - 1
            LLR_left(:,i_left + 1) = f_minsum( LLR(:,i_left + 1), LLR(:,i_left + 1 + N/2) );
        end
        
        [y_hat, ind_ord_left, PM, counter, beta_left] = decode_SSCL(LLR_left, PM, in, counter, y_hat, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array);
        
        LLR = [LLR; LLR];
        LLR = LLR(ind_ord_left(1:l),:);
        
        LLR_right = zeros(l,N/2);
        for i_right = 0:N/2 - 1
            LLR_right(:,i_right + 1) = g_minsum( beta_left(:,i_right + 1), LLR(:,i_right + 1), LLR(:,i_right + 1 + N/2) );
        end
        [y_hat, ind_ord_right, PM, counter, beta_right] = decode_SSCL(LLR_right, PM, in+N/2, counter, y_hat, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array);
        
        ind_ord = [ind_ord_left(1:l) ind_ord_left(1:l)];
        ind_ord = ind_ord(ind_ord_right);
        
        % Updating beta
        beta_left_temp = [beta_left; beta_left];
        beta_left = beta_left_temp(ind_ord_right(1:l),:);
               
        beta = zeros(l,N);
        for i_beta = 0:N/2 - 1
            beta(:,i_beta + 1) = mod( beta_left(:,i_beta + 1) + beta_right(:,i_beta + 1), 2 );
        end
        for i_beta = N/2:N - 1
            beta(:,i_beta + 1) = beta_right(:,i_beta + 1 - N/2);
        end
        
    end
end