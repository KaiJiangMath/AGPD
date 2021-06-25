function u = obtainIniVals(phase, ncpt, dim)
	
	DimPhy = dim(1);
	DimCpt = dim(2);
    u = zeros(ncpt);
    initval = markGlobalPattern(phase);
	kindex = initval(:, 1:DimCpt);
	[nr, nc] = size(kindex);
	for i=1:1:nr
		for j=1:1:nc
			if (kindex(i,j)<0) 
				kindex(i,j)=kindex(i,j)+ncpt(j);
			end
		end
	end
	kindex = kindex+1;
    if DimCpt==4
		for k=1:1:nr
		   u(kindex(k,1),kindex(k,2),kindex(k,3),kindex(k,4)) = initval(k, DimCpt+1);
		end
	elseif DimCpt==6
		for k=1:1:nr
			u(kindex(k,1),kindex(k,2),kindex(k,3),kindex(k,4),kindex(k,5),kindex(k,6)) = initval(k, DimCpt+1);
		end
	elseif DimCpt==3
		for k=1:1:nr
		   u(kindex(k,1),kindex(k,2),kindex(k,3)) = initval(k, DimCpt+1);
		end
	elseif DimCpt==2
		for k=1:1:nr
		   u(kindex(k,1),kindex(k,2)) = initval(k, DimCpt+1);
		end
	elseif DimCpt==1
		u = zeros(1,ncpt(1));
		for k=1:1:nr
		   u(kindex(k)) = initval(k, DimCpt+1);
		end
    end
    
end
