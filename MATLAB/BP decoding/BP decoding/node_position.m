function node_position(len, in, data_indices, data_indices_sorted, frozen_indices, frozen_indices_sorted)

	global Rate_0_node_pos_array;
	global Rate_1_node_pos_array;
	global SPC_node_pos_array;
	global Rep_node_pos_array;
	global Rep_SPC_node_pos_array;

	ind_vec = in:in+len-1;
    
	if ( is_vec_member(ind_vec, frozen_indices_sorted) )
		Rate_0_node_pos_array(in, length(ind_vec)) = Rate_0_node_pos_array(in, length(ind_vec)) + 1;

	elseif ( is_vec_member(ind_vec, data_indices_sorted) )
		Rate_1_node_pos_array(in, length(ind_vec)) = Rate_1_node_pos_array(in, length(ind_vec)) + 1;

	elseif ( ismember(in, frozen_indices_sorted) && is_vec_member(ind_vec(2:end), data_indices_sorted) )
		SPC_node_pos_array(in, length(ind_vec)) = SPC_node_pos_array(in, length(ind_vec)) + 1;

	elseif ( ismember(ind_vec(end:end), data_indices_sorted) && is_vec_member(ind_vec(1:end-1),frozen_indices_sorted) )
		Rep_node_pos_array(in, length(ind_vec)) = Rep_node_pos_array(in, length(ind_vec)) + 1;

	elseif ( is_vec_member(ind_vec(1:length(ind_vec)/2-1),frozen_indices_sorted) && ismember(ind_vec(length(ind_vec)/2),data_indices_sorted) &&...
        ismember(ind_vec(length(ind_vec)/2 + 1), frozen_indices_sorted) && is_vec_member(ind_vec((length(ind_vec)/2 + 2):end),data_indices_sorted) )
		Rep_SPC_node_pos_array(in, length(ind_vec)) = Rep_SPC_node_pos_array(in, length(ind_vec)) + 1;
    end
    
    if (len>1)
        node_position(len/2, in, data_indices, data_indices_sorted, frozen_indices, frozen_indices_sorted);
        node_position(len/2, in+len/2, data_indices, data_indices_sorted, frozen_indices, frozen_indices_sorted);
    end
    
end