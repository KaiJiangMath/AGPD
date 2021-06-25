%% find the index about tau and gamma from a hamilton matrix (rectangular);
function ind = find_ind(tau, gamma, ham_matrix)

	ind = 0;
	[nr, nc] = size(ham_matrix);
	for i = 1:1:nr
		x = ham_matrix(i,1);
		y = ham_matrix(i,2);
		if ( abs(tau-x)<1.0e-8 && abs(gamma-y)<1.0e-8 )
			ind = i;
		end
	end
end
