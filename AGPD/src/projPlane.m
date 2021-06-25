function kspace = projPlane(PATTERN, uCplx, rcpBox, ncpt, dim)
   
   dof = prod(ncpt);
   index = 0;

   DimPhy = dim(1);
   DimCpt = dim(2);
   projmat = getprojmat(PATTERN, DimPhy, DimCpt);
   projmat = projmat*rcpBox;

   kspace = zeros(dof, DimPhy+1);

   if DimCpt == 4 && DimPhy == 2
	   for j1=1:1:ncpt(1)
		   if (j1>ncpt(1)/2) k1=j1-ncpt(1);
		   else k1 = j1;
		   end
		   for j2=1:1:ncpt(2)
			   if (j2>ncpt(2)/2) k2=j2-ncpt(2);
			   else k2 = j2;
			   end
			   for j3=1:1:ncpt(3)
				   if (j3>ncpt(3)/2) k3=j3-ncpt(3);
				   else k3 = j3;
				   end
				   for j4=1:1:ncpt(4)
					   if (j4>ncpt(4)/2) k4=j4-ncpt(4);
					   else k4 = j4;
					   end
					   index = index+1;
					   kk = [k1, k2, k3, k4];
					   kk = kk-1;
					   klow = projmat*(kk');
					   kspace(index,1)=klow(1);
					   kspace(index,2)=klow(2);
					   kspace(index,3)=uCplx(j1,j2,j3,j4);
				   end
			   end
		   end
	   end
   elseif DimCpt == 3 && DimPhy == 2
	   for j1=1:1:ncpt(1)
		   if (j1>ncpt(1)/2) k1=j1-ncpt(1);
		   else k1 = j1;
		   end
		   for j2=1:1:ncpt(2)
			   if (j2>ncpt(2)/2) k2=j2-ncpt(2);
			   else k2 = j2;
			   end
			   for j3=1:1:ncpt(3)
				   if (j3>ncpt(3)/2) k3=j3-ncpt(3);
				   else k3 = j3;
				   end
				   index = index+1;
				   kk = [k1, k2, k3];
				   kk = kk-1;
				   klow = projmat*(kk');
				   kspace(index,1)=klow(1);
				   kspace(index,2)=klow(2);
				   kspace(index,3)=uCplx(j1,j2,j3);
			   end
		   end
	   end

   elseif DimCpt == 2 && DimPhy == 2
	   for j1=1:1:ncpt(1)
		   if (j1>ncpt(1)/2) k1=j1-ncpt(1);
		   else k1 = j1;
		   end
		   for j2=1:1:ncpt(2)
			   if (j2>ncpt(2)/2) k2=j2-ncpt(2);
			   else k2 = j2;
			   end
			   index = index+1;
			   kk = [k1, k2];
			   kk = kk-1;
			   klow = projmat*(kk');
			   kspace(index,1)=klow(1);
			   kspace(index,2)=klow(2);
			   kspace(index,3)=uCplx(j1,j2);
		   end
	   end
   end

end
