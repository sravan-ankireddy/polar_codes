% function to check whether the vector belongs to the set
function res = is_member(ind, indices_sorted)

	res = 0;

	N = length(ind);

	% Staring index
	in = ind(1);

	index = 1:length(indices_sorted);
	
	% Position of starting index
	[~,pos] = is_member(in,indices_sorted);

	if (pos > 0)

		% beginning index
		i1 = index(pos);
		% ending index
		i2 = i1 + N -1;

		if ( i2 <= length(indices_sorted) )

			if ( indices_sorted(i2) == in + N -1 )

				res = 1;

			end

		end
	end
end