
function projmat = getprojmat(PATTERN, DimPhy, DimCpt)

	projmat = zeros(DimPhy, DimCpt);
	if DimCpt == DimPhy  %%%% for crystals
		for j=1:1:DimCpt projmat(j,j)=1.0; end
	elseif DimCpt > DimPhy  %%%% for quasicrystals
		if strcmp(PATTERN,'12fold') || strcmp(PATTERN, 'LQ6') || strcmp(PATTERN, 'LQS6') ||...
			strcmp(PATTERN, 'C3') || strcmp(PATTERN, '12i6o') || strcmp(PATTERN, '8i10o') ||...
		   	strcmp(PATTERN, 'sq') || strcmp(PATTERN, 'squ') || strcmp(PATTERN, 'sqv') ||...
			strcmp(PATTERN, 'sqw') || strcmp(PATTERN, 'Ls')
			nfold = 12;
			for j=1:1:DimCpt
				projmat(1,j)=cos(2*pi*(j-1)/nfold);
				projmat(2,j)=sin(2*pi*(j-1)/nfold);
			end
		elseif strcmp(PATTERN,'10fold')
			nfold = 10;
			for j=1:1:DimCpt
				projmat(1,j)=cos(2*pi*(j-1)/nfold);
				projmat(2,j)=sin(2*pi*(j-1)/nfold);
			end
		elseif strcmp(PATTERN,'icosahedron')
			dimcpt = 6;
			dimphy = 3;
			t = 2*cos(pi/5);
			v1 = [1, 0, 0];
			v2 = 0.5*[t, 1, t-1];
			v4 = 0.5*[t,-1, 1-t];
			v7 = 0.5*[1,t-1,-t];
			v14 = [0, 1, 0];
			v15 = [0, 0, 1];
			projmat = [v1; v2; v4; v7; v14; v15 ];
			projmat = projmat';
		end
	end
end
