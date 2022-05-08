clc
clear

f_c = 1*10^9;
f_s = 5*10^9; %Typical maximum cpu clock
fraction_list = linspace(0.00001, 1, 300);
weights = 38469 * 32;
data_size = 3925000 * 8;
Budget = 70*10^3;
cpu_parameter = 2*10^(-28);
pay_off = 10^(-9);
privacy_list = [7500, 7600, 7700, 7800, 7900, 8000, 8100, 8200, 8300, 8400];
gamma = 0.001;
NBI_resolution = 300;
PENALTY_INCREASE = 10;
threshold = 0.001;
K = 5;
N = 10;
rho = 1000;
fraction_outcomes = zeros(N, 1);


%%%%%%NBS parameter%%%%%%%
Disagree_point_1 = 0;
Disagree_point_2 = 0;
CONSTANT_small = 1;
CONSTANT_big = 999999999;
threshold_distant = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% NBI method %%%%%%%%
privacy_count = 1;
U_1_baseline = zeros(N, 1);
U_1_0p1 = zeros(N,1);
U_1_Op5 = zeros(N, 1);
U_1_Op7 = zeros(N, 1);


for privacy_coefficient = privacy_list
    U_1_ideal = privacy_coefficient/(log(2) * (weights * data_size * cpu_parameter * f_c^2 - pay_off*f_c)) - 1;
    U_2_ideal = 0;
    U_1 = @(alpha) pay_off*f_c - alpha * weights * data_size * cpu_parameter * f_c^2 ...
        +privacy_coefficient * log2(1 + alpha);
    U_2 = @(alpha) Budget - gamma * (1 - alpha) * weights * data_size * cpu_parameter * f_s^2 ...
        -(1 - gamma)* (alpha * weights *K* data_size/f_c + (1 - alpha)*K * weights * data_size/f_s + rho* log2(1+K/N));
    U_1_baseline(privacy_count) = U_1(1);
    U_1_0p1(privacy_count) = U_1(0.1);
    U_1_0p5(privacy_count) = U_1(0.5);
    U_1_0p7(privacy_count) = U_1(0.7);
    
    
    [U_1_Pareto, U_2_Pareto, Pareto_optimal_points] = NBI_4(U_1, U_2, NBI_resolution, PENALTY_INCREASE, threshold, U_1_ideal, U_2_ideal);
    
    %%%%%%NBS RUN%%%%%%%%
    [NBS_client, NBS_server, tangential, CONSTANT, NBS_INPUT, NBS_OUTPUT] = bisection(U_1_Pareto,CONSTANT_small, CONSTANT_big, U_2_Pareto, Disagree_point_1, Disagree_point_2, threshold_distant);
    % [NBS_client, NBS_server, tangential, CONSTANT, NBS_INPUT, NBS_OUTPUT] = bisection(Utility_client,CONSTANT_small, CONSTANT_big, Utility_server, Disagree_point_1, Disagree_point_2, threshold_distant);
    fraction_outcomes(privacy_count) = fraction_list(tangential+1);
    privacy_count = privacy_count + 1;
end

avg_fraction = mean(fraction_outcomes);
avg_fraction = 0.3; %rounding

U_1_proposed = pay_off*f_c - avg_fraction * weights * data_size * cpu_parameter * f_c^2 ...
    +privacy_list.* log2(1 + avg_fraction);
U_2_0p1 = U_2(0.1);
U_2_baseline = U_2(1.0);

U_2_proposed = Budget - gamma * (1 - avg_fraction) * weights * data_size * cpu_parameter * f_s^2 ...
    -(1 - gamma)* (avg_fraction * weights *K* data_size/f_c + (1 - avg_fraction)*K * weights * data_size/f_s + rho* log2(1+K/N));

ref = [1 2 3];

U = [sum(U_1_0p1), sum(U_1_proposed), sum(U_1_baseline)];
b = bar(ref, U, 'FaceColor','flat');
% bar(F, U, 'k')
grid on;
% xlabel('\alpha')
ylabel('Sum of the utilities of clients')
name = {'\alpha = 0.1', '\alpha = 0.3', '\alpha = 1.0'};
set(gca, 'xticklabel', name)


