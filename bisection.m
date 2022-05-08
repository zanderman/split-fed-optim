function [NBS_client, NBS_server, tangential, CONSTANT, NBS_INPUT, PARA_OUTPUT]  = bisection(Utility_client,CONSTANT_small, CONSTANT_big, Utility_server, Disagree_point_1, Disagree_point_2, threshold_distant)

L = length(Utility_client);
Utility_client = Utility_client(2 : L-1);
Utility_server = Utility_server(2 : L-1); %%Preprocessing
distance = zeros(1, length(Utility_client));
CONSTANT = (CONSTANT_small + CONSTANT_big)/2;
PARA_OUTPUT = zeros(1, length(distance));

while true

    para = @(x) (CONSTANT + Disagree_point_1 * x - Disagree_point_1 * Disagree_point_2) ...
        /(x - Disagree_point_2); %%parameterized parabola 
   
   
%     NASH_OUTPUT = para(Utility_client); %the output of the parabola given NASH_INPUT, which is the energy consumption
    PARA_OUTPUT = CONSTANT./Utility_client;
    
    distance = PARA_OUTPUT - Utility_server; %distance of y-axis between #of global iterations and NASH_OUTPUT
    
    if min(distance) > threshold_distant %the parabola is ahead of the Pareto boudary and it is not close enough 
        CONSTANT_big = CONSTANT;
        CONSTANT = (CONSTANT_small + CONSTANT)/2;
        
    elseif min(distance) < 0 %the parabola is below the Pareto boundary. The constant is too small 
        CONSTANT_small = CONSTANT;
        CONSTANT = (CONSTANT_big + CONSTANT)/2;
        
    elseif min(distance) < threshold_distant && min(distance) > 0 %the parabola is below the Pareto boundary and it is close enough 
        [foo, tangential] = min(distance);
%         tangential = round(tangential);
        NBS_client = Utility_client(tangential);
        NBS_server = Utility_server(tangential);
        break
    end
end
NBS_INPUT = Utility_client;
end