function uc = mymap(uc_cs, uc_rf, mapflag)
	
	n_cs = size(uc_cs);
	n_rf = size(uc_rf);

	tmp_rf = fftshift(uc_rf);
	tmp_cs = fftshift(uc_cs);

	if length(n_cs) == 4
		if	strcmp('cs2rf', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			k2 = n_rf(2)/2-n_cs(2)/2;
			k3 = n_rf(3)/2-n_cs(3)/2;
			k4 = n_rf(4)/2-n_cs(4)/2;

			tmp_rf(k1+1:1:k1+n_cs(1), k2+1:1:k2+n_cs(2), k3+1:1:k3+n_cs(3), k4+1:1:k4+n_cs(4)) = tmp_cs; 
			uc = fftshift(tmp_rf);

		elseif strcmp('rf2cs', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			k2 = n_rf(2)/2-n_cs(2)/2;
			k3 = n_rf(3)/2-n_cs(3)/2;
			k4 = n_rf(4)/2-n_cs(4)/2;
			
			tmp_cs = tmp_rf(k1+1:1:k1+n_cs(1), k2+1:1:k2+n_cs(2), k3+1:1:k3+n_cs(3), k4+1:1:k4+n_cs(4));
			uc = fftshift(tmp_cs);
		end
	end

	if length(n_cs) == 3
		if	strcmp('cs2rf', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			k2 = n_rf(2)/2-n_cs(2)/2;
			k3 = n_rf(3)/2-n_cs(3)/2;

			tmp_rf(k1+1:1:k1+n_cs(1), k2+1:1:k2+n_cs(2), k3+1:1:k3+n_cs(3)) = tmp_cs; 
			uc = fftshift(tmp_rf);

		elseif strcmp('rf2cs', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			k2 = n_rf(2)/2-n_cs(2)/2;
			k3 = n_rf(3)/2-n_cs(3)/2;
			
			tmp_cs = tmp_rf(k1+1:1:k1+n_cs(1), k2+1:1:k2+n_cs(2), k3+1:1:k3+n_cs(3));
			uc = fftshift(tmp_cs);
		end
	end

	if length(n_cs) == 2
		if	strcmp('cs2rf', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			k2 = n_rf(2)/2-n_cs(2)/2;

			tmp_rf(k1+1:1:k1+n_cs(1), k2+1:1:k2+n_cs(2)) = tmp_cs; 

			uc = fftshift(tmp_rf);

		elseif strcmp('rf2cs', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			k2 = n_rf(2)/2-n_cs(2)/2;
			
			tmp_cs = tmp_rf(k1+1:1:k1+n_cs(1), k2+1:1:k2+n_cs(2));
			uc = fftshift(tmp_cs);
		end
	end

	if length(n_cs) == 1
		if	strcmp('cs2rf', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			tmp_rf(k1+1:1:k1+n_cs(1)) = tmp_cs; 
			uc = fftshift(tmp_rf);

		elseif strcmp('rf2cs', mapflag)
			k1 = n_rf(1)/2-n_cs(1)/2;
			tmp_cs = tmp_rf(k1+1:1:k1+n_cs(1));
			uc = fftshift(tmp_cs);
		end
	end

end
