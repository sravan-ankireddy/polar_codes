function [LLR, codeword, beta] = decode(LLR, in, ind_lev, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array, Rep_SPC_node_pos_array)
	
	% No of levels
	n0 = size(LLR,1);
    
	% Length of current sub-vector
	N = 2^(n0 - ind_lev);
    
    codeword = zeros(1,N);

	% LLR corresponding to current level
	LLR_lev = LLR(ind_lev,in:in+N-1);
%     disp(LLR_lev);

	% Sub-indices vector
	ind_vec = in:in+N-1;

	if ( Rate_0_node_pos_array(in, length(ind_vec)) > 0 )
		codeword(LLR_lev > 0) = 0;
		codeword(LLR_lev < 0) = 0;
		beta = codeword;

	elseif ( Rate_1_node_pos_array(in, length(ind_vec)) > 0 )
		codeword(LLR_lev > 0) = 0;
		codeword(LLR_lev < 0) = 1;
		beta = codeword;
		codeword = encode(codeword);

	elseif ( SPC_node_pos_array(in, length(ind_vec)) > 0 )
		codeword(LLR_lev > 0) = 0;
		codeword(LLR_lev < 0) = 1;
		parity_bit = mod(sum(codeword), 2);

		% position of least reliable bit
		[~, pos] = min(abs(LLR_lev));

		% Parrity Check correction
		codeword(pos) = mod( (codeword(pos) + parity_bit), 2 );
		beta = codeword;
		codeword = encode(codeword);

	elseif ( Rep_node_pos_array(in, length(ind_vec)) > 0 )
		if (sum(LLR_lev) > 0)
			codeword(LLR_lev > 0) = 0;
			codeword(LLR_lev < 0) = 0;
		else
			codeword(LLR_lev > 0) = 1;
			codeword(LLR_lev < 0) = 1;
		end

		beta = codeword;
		codeword = encode(codeword);

	else

		% Left fork
		LLR_lev_left = zeros(1,N/2);
		for i_left = 0:N/2 - 1
			LLR_lev_left(i_left + 1) = f_minsum( LLR_lev(i_left + 1), LLR_lev(i_left + 1 + N/2) );
        end

        LLR(ind_lev + 1, in:in + N/2 - 1) = LLR_lev_left;

		[LLR, codeword_left, beta_left] = decode(LLR, in, ind_lev + 1, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array, Rep_SPC_node_pos_array);

		% Right fork
		LLR_lev_right = zeros(1,N/2);
		for i_right = 0:N/2 - 1
			LLR_lev_right(i_right + 1) = g_minsum( beta_left(i_right + 1), LLR_lev(i_right + 1), LLR_lev(i_right + 1 + N/2) );
        end

		LLR(ind_lev + 1, in + N/2 :in + N - 1) = LLR_lev_right;

		[LLR, codeword_right, beta_right] = decode(LLR, in + N/2, ind_lev + 1, Rate_0_node_pos_array, Rate_1_node_pos_array, SPC_node_pos_array, Rep_node_pos_array, Rep_SPC_node_pos_array);

		codeword = [codeword_left codeword_right];

		% Updating the beta information
        beta = zeros(1,N);
		for i_beta = 0:N/2 - 1
			beta(i_beta + 1) = mod( beta_left(i_beta + 1) + beta_right(i_beta + 1), 2 );
		end
		for i_beta = N/2:N - 1
			beta(i_beta + 1) = beta_right(i_beta + 1 - N/2);
		end
    end
    
    
end