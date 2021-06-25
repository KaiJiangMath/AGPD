function potenCplx = leibler_potential(q2)

	epsilon = 1;
	q2(1) = 1;
	f = 0.5;

	switch epsilon
	case 1
		min_q = 1.945608;
		min_poten = 20.989737; 
		min_dd_poten = 14.570483;
	case 2
		min_q = 1.213181;
		min_poten = 21.792811;
		min_dd_poten = 37.671549;
	case 3
		min_q = 0.846117;
		min_poten = 22.465837;
		min_dd_poten = 78.667164;
	case 4
		min_q = 0.643937; 
		min_poten = 22.814863;
		min_dd_poten = 137.454230;
	case 5
		min_q = 0.518450;
		min_poten = 23.004400;
		min_dd_poten = 213.642168;
	case 6
		min_q = 0.433458;
		min_poten = 23.115957;
		min_dd_poten = 307.103783;
	case 7
		min_q = 0.372264;
		min_poten = 23.186381;
		min_dd_poten = 417.662503;
	case 8
		min_q = 0.326069;
		min_poten = 23.233420;
		min_dd_poten = 545.656519;
	case 9
		min_q = 0.290072;
		min_poten = 23.266292;
		min_dd_poten = 690.553261;
	case 10
		min_q = 0.261275;
		min_poten = 23.290122;
		min_dd_poten = 851.941738;
	case 11
		min_q = 0.237577;
		min_poten = 23.307926;
		min_dd_poten = 1031.436730;
	case 12
		min_q = 0.217879;
		min_poten = 23.321566;
		min_dd_poten = 1227.010452;
	case 13
		min_q = 0.201181;
		min_poten = 23.332241;
		min_dd_poten = 1439.815922;
	case 14
		min_q = 0.186882;
		min_poten = 23.340751;
		min_dd_poten = 1668.954307;
	case 15
		min_q = 0.174384;
		min_poten = 23.347638;
		min_dd_poten = 1918.110086;
	case 16
		min_q = 0.163485;
		min_poten = 23.353293;
		min_dd_poten = 2183.242456;
	case 17
		min_q = 0.153886;
		min_poten = 23.357990;
		min_dd_poten = 2464.658558;
	case 18
		min_q = 0.145387;
		min_poten = 23.361931;
		min_dd_poten = 2761.106938;
	case 19
		min_q = 0.137787;
		min_poten = 23.365277;
		min_dd_poten = 3073.684158;
	case 20
		min_q = 0.130888;
		min_poten = 23.368133;
		min_dd_poten = 3407.195159;
	end
	
	eps2 = epsilon.^2;	eps4 = eps2.^2;
	q2 = q2 * (min_q^2);
	q4 = q2.^2;

	S11 = 2./(eps4.*q4) .* ( f.*eps2.*q2 + exp(-f.*eps2.*q2) - 1 );
	S22 = 2./q4 .* ( (1-f).*q2 + exp(-(1-f).*q2) - 1 );
	S12 = 1./(eps2.*q4) .* (1 - exp(-f.*eps2.*q2)) .* (1 - exp(-(1-f).*q2));

	potenCplx = (S11 + S22 + 2*S12) ./ (S11.*S22 - S12.*S12);
	potenCplx = 8/min_dd_poten * (potenCplx - min_poten);
	potenCplx(1) = 0;
end
