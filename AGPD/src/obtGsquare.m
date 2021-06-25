%% generate the Fourier k with the given computational box
function Gsquare = obtGsquare(ncpt, rcpBox, dim, projmat)

	%% prepare
    format long;
	DimPhy = dim(1);
	DimCpt = dim(2);
	Gsquare = zeros(ncpt);

	%% taking the product
	projBox = projmat * rcpBox;
	[nr, nc] = size(projBox);

	%% generate the Fourier PBk^2
	if ( nc == 2 )
		kk1 = repmat( [0:1:ncpt(1)/2, -ncpt(1)/2+1:1:-1].', [1,ncpt(2)] );
		kk2 = repmat( [0:1:ncpt(2)/2, -ncpt(2)/2+1:1:-1], [ncpt(1),1] );
		%% compute the Fourier (PBk)^2
		for j = 1:1:nr
			temp = projBox(j,1) * kk1 + projBox(j,2) * kk2;
			Gsquare = Gsquare + temp.^2;
		end
	elseif ( nc == 3 )
		kk1 = repmat( [0:1:ncpt(1)/2, -ncpt(1)/2+1:1:-1].', [1,ncpt(2),ncpt(3)] );
		kk2 = repmat( [0:1:ncpt(2)/2, -ncpt(2)/2+1:1:-1], [ncpt(1),1,ncpt(3)] );
		kk3 = repmat( [0:1:ncpt(3)/2, -ncpt(3)/2+1:1:-1].', [1,ncpt(1),ncpt(2)] );
		kk3 = permute( kk3, [2,3,1] );
		%% compute the Fourier (PBk)^2
		for j = 1:1:nr
			temp = projBox(j,1) * kk1 + projBox(j,2) * kk2 + projBox(j,3) * kk3;
			Gsquare = Gsquare + temp.^2;
		end
	elseif ( nc == 4 )
		kk1 = repmat( [0:1:ncpt(1)/2, -ncpt(1)/2+1:1:-1].', [1,ncpt(2),ncpt(3),ncpt(4)] );
		kk2 = repmat( [0:1:ncpt(2)/2, -ncpt(2)/2+1:1:-1], [ncpt(1),1,ncpt(3),ncpt(4)] );
		kk3 = repmat( [0:1:ncpt(3)/2, -ncpt(3)/2+1:1:-1].', [1,ncpt(1),ncpt(2),ncpt(4)] );
		kk3 = permute( kk3, [2,3,1,4] );
		kk4 = repmat( [0:1:ncpt(4)/2, -ncpt(4)/2+1:1:-1].', [1,ncpt(1),ncpt(2),ncpt(3)] );
		kk4 = permute( kk4, [2,3,4,1] );
		%% compute the Fourier (PBk)^2
		for j = 1:1:nr
			temp = projBox(j,1) * kk1 + projBox(j,2) * kk2 + projBox(j,3) * kk3 + projBox(j,4) * kk4;
			Gsquare = Gsquare + temp.^2;
		end
	end
end
