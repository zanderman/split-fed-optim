function [U_1_Pareto, U_2_Pareto, Pareto_optimal_points] = NBI_4(U_1, U_2, NBI_resolution, PENALTY_INCREASE, threshold, U_1_ideal, U_2_ideal)
options = optimoptions(@fmincon, 'Display', 'none');
lb = 0;
ub = 1;
alpha_0 = 0.3;

NBI_array =linspace(0, 1, NBI_resolution);
NBI_increase = NBI_array(2) - NBI_array(1);
NBI = 0; %beta_1

U_1_Pareto = zeros(1, length(NBI_array));
U_2_Pareto = zeros(1, length(NBI_array));
Pareto_optimal_points = zeros(length(NBI_array), 1);

index = 1;

while index < NBI_resolution+1
    
    PENALTY = 1;
    count = 1;
    
    while true
        if count == 1
            fnc = @(alpha) -(U_2(alpha) - U_2(U_2_ideal))/(U_2(U_1_ideal) -U_2(U_2_ideal)) + NBI + PENALTY * ...
                (1/2 - 3/2*NBI + (U_2(alpha) - U_2(U_2_ideal))/(U_2(U_1_ideal) -U_2(U_2_ideal)) ...
                - 1/2 * (U_1(alpha) - U_1(U_1_ideal))/(U_1(U_2_ideal) - U_1(U_1_ideal)) )^2;
            Present_point = fmincon(fnc, alpha_0, [], [], [] ,[], lb, ub, [], options);
        else
            Present_point = Next_point;
        end
        PENALTY = PENALTY_INCREASE * PENALTY;
            fnc = @(alpha) -(U_2(alpha) - U_2(U_2_ideal))/(U_2(U_1_ideal) -U_2(U_2_ideal)) + NBI + PENALTY * ...
                (1/2 - 3/2*NBI + (U_2(alpha) - U_2(U_2_ideal))/(U_2(U_1_ideal) -U_2(U_2_ideal)) ...
                - 1/2 * (U_1(alpha) - U_1(U_1_ideal))/(U_1(U_2_ideal) - U_1(U_1_ideal)) )^2;
        Next_point = fmincon(fnc, Present_point, [], [], [] ,[], lb, ub, [], options); %start from the previous optimal, it is much faster
        koo = Present_point - Next_point;
        if abs(koo) < threshold
            Pareto_optimum = Next_point;
            U_1_Pareto(index) = U_1(Pareto_optimum);
            U_2_Pareto(index) = U_2(Pareto_optimum);
            Pareto_optimal_points(index) = Pareto_optimum;
            break
        end
        count = count + 1;
    end
    NBI = NBI + NBI_increase;
    index = index + 1;
end