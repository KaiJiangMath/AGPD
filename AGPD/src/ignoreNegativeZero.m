%% ignore the negative zero
%% only consider the positive zero
function val = ignoreNegativeZero(val)
	if ( abs(val) == val )
		val = abs(val);
	end
end
