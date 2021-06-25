%% screen these candidate phases
function data = deal(data)
	ind = find(data(:,4)~=1);
	data(ind,3) = 0;
end
