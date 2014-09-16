%% Antenna and propagation mm1 

% Excercise 2.4 Find the HPBW and FNBW  in radians and degress for the
% following normalized radiation intensities 
% symbolic link 
syms theta;

% a U(theta) = cos(theta)

ans_a_hpbw = 2*solve(cos(theta) == 0.5);
ans_a_fnbw = 2*solve(cos(theta) == 0);
ans_a_hpbw_degree = int16(radtodeg(ans_a_hpbw));
ans_a_fnbw_degree = int16(radtodeg(ans_a_fnbw));

disp('HPBW for ex a')
disp(ans_a_hpbw(1,1))
disp(ans_a_hpbw_degree(1,1))
disp('FNBW for ex a')
disp(ans_a_fnbw)
disp(ans_a_fnbw_degree)
% b U(theta) = cos^2(theta)
ans_b_hpbw = 2*solve(cos(theta)^2 == 0.5);
ans_b_fnbw = 2*solve(cos(theta)^2 == 0);
ans_b_hpbw_degree = int16(radtodeg(ans_b_hpbw));
ans_b_fnbw_degree = int16(radtodeg(ans_b_fnbw));

disp('HPBW for ex b')
disp(ans_b_hpbw(1,1))
disp(ans_b_hpbw_degree(1,1))
disp('FNBW for ex b')
disp(ans_b_fnbw)
disp(ans_b_fnbw_degree)

% c U(theta) = cos(2*theta)
ans_c_hpbw = 2*solve(cos(2*theta) == 0.5);
ans_c_fnbw = 2*solve(cos(2*theta) == 0);
ans_c_hpbw_degree = int16(radtodeg(ans_c_hpbw));
ans_c_fnbw_degree = int16(radtodeg(ans_c_fnbw));
disp('HPBW for ex c')
disp(ans_c_hpbw(1,1))
disp(ans_c_hpbw_degree(1,1))
disp('FNBW for ex c')
disp(ans_c_fnbw)
disp(ans_c_fnbw_degree)

% d U(theta) = cos^2(2*theta)
ans_d_hpbw = 2*solve(cos(2*theta)^2 == 0.5);
ans_d_fnbw = 2*solve(cos(2*theta)^2 == 0);
ans_d_hpbw_degree = int16(radtodeg(ans_d_hpbw));
ans_d_fnbw_degree = int16(radtodeg(ans_d_fnbw));
disp('HPBW for ex d')
disp(ans_d_hpbw(1,1))
disp(ans_d_hpbw_degree(1,1))
disp('FNBW for ex d')
disp(ans_d_fnbw)
disp(ans_d_fnbw_degree)
% e U(theta) = cos(3*theta)
ans_e_hpbw = 2*solve(cos(3*theta) == 0.5);
ans_e_fnbw = 2*solve(cos(3*theta) == 0);
ans_e_hpbw_degree = int16(radtodeg(ans_e_hpbw));
ans_e_fnbw_degree = int16(radtodeg(ans_e_fnbw));
disp('HPBW for ex e')
disp(ans_e_hpbw(1,1))
disp(ans_e_hpbw_degree(1,1))
disp('FNBW for ex e')
disp(ans_e_fnbw)
disp(ans_e_fnbw_degree)
% f U(theta) = cos^2(3*theta)
ans_f_hpbw = 2*solve(cos(3*theta)^2 == 0.5);
ans_f_fnbw = 2*solve(cos(3*theta)^2 == 0);
ans_f_hpbw_degree = int16(radtodeg(ans_f_hpbw));
ans_f_fnbw_degree = int16(radtodeg(ans_f_fnbw));
disp('HPBW for ex f')
disp(ans_f_hpbw(1,1))
disp(ans_f_hpbw_degree(1,1))
disp('FNBW for ex f')
disp(ans_f_fnbw)
disp(ans_f_fnbw_degree)