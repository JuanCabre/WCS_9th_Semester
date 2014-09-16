%% Antenna and propagation mm1 

% Excercise 2.4 Find the HPBW and FNBW  in radians and degress for the
% following normalized radiation intensities 
% symbolic link 
syms theta;

% a U(theta) = cos(theta)

ans_a_hpbw = solve(cos(theta) == 0.5);
ans_a_fnbw = solve(cos(theta) == 0);
disp('HPBW for ex a')
disp(ans_a_hpbw)
disp('FNBW for ex a')
disp(ans_a_fnbw)

% b U(theta) = cos^2(theta)
ans_b_hpbw = solve(cos(theta)^2 == 0.5);
ans_b_fnbw = solve(cos(theta)^2 == 0);
disp('HPBW for ex b')
disp(ans_b_hpbw)
disp('FNBW for ex b')
disp(ans_b_fnbw)

% c U(theta) = cos(2*theta)
ans_c_hpbw = solve(cos(2*theta) == 0.5);
ans_c_fnbw = solve(cos(2*theta) == 0);
disp('HPBW for ex c')
disp(ans_c_hpbw)
disp('FNBW for ex c')
disp(ans_c_fnbw)

% d U(theta) = cos^2(2*theta)
ans_d_hpbw = solve(cos(2*theta)^2 == 0.5);
ans_d_fnbw = solve(cos(2*theta)^2 == 0);
disp('HPBW for ex d')
disp(ans_d_hpbw)
disp('FNBW for ex d')
disp(ans_d_fnbw)

% e U(theta) = cos(3*theta)
ans_e_hpbw = solve(cos(3*theta) == 0.5);
ans_e_fnbw = solve(cos(3*theta) == 0);
disp('HPBW for ex e')
disp(ans_e_hpbw)
disp('FNBW for ex e')
disp(ans_e_fnbw)

% f U(theta) = cos^2(3*theta)
ans_f_hpbw = solve(cos(3*theta)^2 == 0.5);
ans_f_fnbw = solve(cos(3*theta)^2 == 0);
disp('HPBW for ex f')
disp(ans_f_hpbw)
disp('FNBW for ex f')
disp(ans_f_fnbw)