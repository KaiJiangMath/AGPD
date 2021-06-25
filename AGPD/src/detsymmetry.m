function sflag = detsymmetry(uc, u0)
%%%  u0: principle symmetric wave vectors
%%%  uc: the results required to determine 

	NdetSymm = sum(uc(:)~=0);
    [us, uind] = sort(abs(uc(:)), 'descend');
    uind = uind(1:NdetSymm);
    [u0s, u0ind] = sort(abs(u0(:)), 'descend');
    u0ind = u0ind(1:NdetSymm);
    udet = intersect(uind, u0ind);
    if length(udet) == NdetSymm 
        sflag = 1;
    else
        sflag = 0;
    end
end
