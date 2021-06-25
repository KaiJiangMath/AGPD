
function potenCplx = obtainPotential(PATTERN, ncpt, rcpBox, dim)
    format long;
	DimPhy = dim(1);
	DimCpt = dim(2);
	projmat = getprojmat(PATTERN, DimPhy, DimCpt);
	potenCplx=zeros(ncpt);
	ss = (projmat')*projmat;

	if DimCpt == 3
		for j1=1:1:ncpt(1)
			if (j1>ncpt(1)/2+1) k1=j1-ncpt(1);
			else k1 = j1;
			end
			for j2=1:1:ncpt(2)
				if (j2>ncpt(2)/2+1) k2=j2-ncpt(2);
				else k2 = j2;
				end
				for j3=1:1:ncpt(3)
					if (j3>ncpt(3)/2+1) k3=j3-ncpt(3);
					else k3 = j3;
					end
					kk = [k1, k2, k3];
					kk = kk-1;
					t2 = kk*(rcpBox')*ss*rcpBox*(kk');
%                    potenCplx(j1,j2,j3) = (1-t2);    %% LB potential
					if kk(1) ~= 0 & kk(2) ~= 0 & kk(2) ~= 0 
						potenCplx(j1,j2,j3) = t2 + 1/t2 -2;    %% OK potential
					end
					potenCplx(1,1,1) = 0;    %% OK potential
				end
			end
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif DimCpt == 2
		for j1=1:1:ncpt(1)
			if (j1>ncpt(1)/2+1) k1=j1-ncpt(1);
			else k1 = j1;
			end
			for j2=1:1:ncpt(2)
				if (j2>ncpt(2)/2+1) k2=j2-ncpt(2);
				else k2 = j2;
				end
				kk = [k1, k2];
				kk = kk-1;
				t2 = kk*(rcpBox')*ss*rcpBox*(kk');
				potenCplx(j1,j2) = (1-t2);
			end
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif DimCpt == 1
		for j1=1:1:ncpt(1)
			if (j1>ncpt(1)/2+1) k1=j1-ncpt(1);
			else k1 = j1;
			end
			kk = k1-1;
			t2 = kk*(rcpBox')*ss*rcpBox*(kk');
			potenCplx(j1) = (1-t2);
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif (DimCpt == 4) && (DimPhy == 2)

		for j1=1:1:ncpt(1)
			if (j1>ncpt(1)/2+1) k1=j1-ncpt(1);
			else k1 = j1;
			end
			for j2=1:1:ncpt(2)
				if (j2>ncpt(2)/2+1) k2=j2-ncpt(2);
				else k2 = j2;
				end
				for j3=1:1:ncpt(3)
					if (j3>ncpt(3)/2+1) k3=j3-ncpt(3);
					else k3 = j3;
					end
					for j4=1:1:ncpt(4)
						if (j4>ncpt(4)/2+1) k4=j4-ncpt(4);
						else k4 = j4;
						end
						kk = [k1, k2, k3, k4];
						kk = kk-1;
						t2 = kk*(rcpBox')*ss*rcpBox*(kk');
						potenCplx(j1,j2,j3,j4) = (1-t2);
					end
				end
			end
		end

	end
%%%%%%%%%%%%%%%%%%%%%%%%%%
end

