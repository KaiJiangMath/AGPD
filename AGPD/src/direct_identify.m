%% identify structures by the Hash value
function [similarity, sflag] = direct_identify(PATTERN, file_obj, file_ref)
	%% file_obj: the objective
	%% file_ref: the reference

	N = 512;		% scaling size for identifying
	% the reference value: <=refN, similar; >refN, non-similar
	if strcmp(PATTERN, 'fcc')
		refN = 15000;
	else
		refN = 12000;
	end
	hash_obj = aHash(file_obj, N);
	hash_ref = aHash(file_ref, N);

	%% the similarity degree
	similarity = sum(abs(hash_obj - hash_ref));
    if ( similarity <= refN )
        sflag = 1;
    else
        sflag = 0;
    end
end
