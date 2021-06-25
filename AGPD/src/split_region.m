%% split the computational region
function taur = split_region(taur_tot, i, tau_num)

	% array
	dtau = taur_tot(3);
	ltau = round((taur_tot(2)-taur_tot(1))/dtau) + 1;

	taur = zeros(1,3);	% [start, end, step];
	taur(3) = dtau;
	if ( tau_num <= 1 )
		taur = taur_tot;
		ds = 1;
	else
		ds = floor(ltau/tau_num);
		excess = round(ltau - ds*tau_num);
		if ( i < excess )
			ds = ds + 1;
		end
		taur(1) = taur_tot(1) + dtau*i*ds;
		taur(2) = taur_tot(1) + dtau*(i+1)*ds - dtau;
		if ( i >= excess )
			taur(1:2) = taur(1:2) + dtau*excess;
		end
	end
	if ( abs(taur(1)-taur(2)) < 1e-8 )
		taur(2) = taur(1);
	end
	fprintf('--> %d (%d) \t %d \t [%.6f, %.6f, %.6f]\n', i, tau_num, ds, taur(1), taur(2), dtau);
end
