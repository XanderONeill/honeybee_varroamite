function dydt = HoneyBee_disperse_and_stay(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params)
threshold = 1;

N = length(y)/5;            % number of hives
N_A = sum((y(2:5:end)+y(4:5:end)+y(5:5:end))>threshold);    % number of alive hives
% N_A = N - sum(y(2:5:end)==0);

% type - vector of length N, detailing what category of bee is in each hive
% if type(i) == 1, regular hive
% if type(i) == 2, tolerant hive
% if type(i) == 3, resistant hive

dydt = zeros(N*5, 1);

%% Set Parameters - fixed regardless of colony/hive type
% BB_p - Bee Brood Parameters
% HB_p - Honeybee Parameters
% VM_p - Varroa mite parameters
% SP_p - Space parameters

a_H = BB_p(1);      omega = BB_p(2);    b_B = BB_p(3);      m = BB_p(4);
q_H = HB_p(1);
b_V = VM_p(1);      alpha_B = VM_p(2);  sigma = VM_p(3);    
nu  = SP_p(1);

if rem(t, 12) <= 3
    a_H = 1*a_H/3;
elseif (rem(t,12) > 6) && (rem(t,12) <= 9)
    a_H = 1*a_H/3;
elseif rem(t,12) > 9
    a_H = 0*a_H/9;
end

%% Parameters that depend on hive type (regular, resistant and tolerant)

% type_params(1,:) = a_V        - entry 1 - parameter for tolerant and regular hives
%                               - entry 2 - parameter for resistant bees
% type_params(2,:) = alpha_H    - entry 1 - parameter for regular and resistant bees
%                               - entry 2 - parameter for tolerant bees

%% Model Equations

for i = 1:N

    if types(i) == 1                     % Resident Regular Bees
        a_V = type_params(1,1);
        alpha_H = type_params(1,2);
      
        % if resident (R) is regular, then visiting honeybee classes are
        % tolerant and resistant
        H_R = y(5*find(types == 1) - 3);
        V_R = y(5*find(types == 1) - 2);
        H_1 = y(5*find(types == 2) - 3);
        V_1 = y(5*find(types == 2) - 2);
        H_2 = y(5*find(types == 3) - 3);
        V_2 = y(5*find(types == 3) - 2);
    elseif types(i) == 2                 % Resident Tolerant Bees
        a_V = type_params(1,1);
        alpha_H = type_params(2,2);

        % if resident (R) is tolerant, then visiting honeybee classes are
        % regular and resistant
        H_R = y(5*find(types == 2) - 3);
        V_R = y(5*find(types == 2) - 2);
        H_1 = y(5*find(types == 1) - 3);
        V_1 = y(5*find(types == 1) - 2);
        H_2 = y(5*find(types == 3) - 3);
        V_2 = y(5*find(types == 3) - 2);
    elseif types(i) == 3                 % Resident Resistant Bees
        a_V = type_params(2,1);
        alpha_H = type_params(1,2);

        % if resident (R) is resistant, then visiting honeybee classes are
        % regular and tolerant
        H_R = y(5*find(types == 3) - 3);
        V_R = y(5*find(types == 3) - 2);
        H_1 = y(5*find(types == 1) - 3);
        V_1 = y(5*find(types == 1) - 2);
        H_2 = y(5*find(types == 2) - 3);
        V_2 = y(5*find(types == 2) - 2);
    end

    H_T = y(5*i-3) + y(5*i-1) + y(5*i-0);   %total number of honeybees in hive i
    V_T = y(5*i-2);                         %total varroa in hive i

    Burden = V_T./max(1, H_T);
    Burden(isinf(Burden)) = 0;
    Burden(isnan(Burden)) = 0;

    if H_T > 10000
        b_H = 0;
    else
        b_H = HB_p(2);
    end

    % RESIDENT BROOD
    % Brood eggs are laid by the 1 resident honeybee queen, but any
    % honeybees in the hive can contribute towards reproduction
    dydt(5*i-4) =   a_H*H_T^2/(omega^2 + H_T^2)...
                    - b_B*y(5*i-4)...
                    - m*y(5*i-4)...
                    - alpha_B*V_T*y(5*i-4)/(sigma + y(5*i-4));

    % RESIDENT HONEYBEES
    % We assume all brood progress to a honeybee of the resident class
    if H_T > threshold
        dydt(5*i-3) =   m*y(5*i-4)...
                        - q_H*H_T*y(5*i-3)...
                        - alpha_H*Burden*y(5*i-3)...
                        - b_H*y(5*i-3)...
                        - nu*y(5*i-3)...
                        + nu/(N_A)*sum(H_R);
    else
        dydt(5*i-3) =   m*y(5*i-4)...
                        - q_H*H_T*y(5*i-3)...
                        - alpha_H*Burden*y(5*i-3)...
                        - b_H*y(5*i-3);
    end

    % RESIDENT VARROA MITES
    % We assume that if there are no honeybees of the resident type in the
    % resident hive, then no other bees would want to visit, nor would any
    % varroa mites be able to leave (as there are no resident bees to carry
    % them). Hence, all movement terms are multiplied by min(y(5*i-3), 1)
    if H_T > threshold
        dydt(5*i-2) =   a_V*V_T*y(5*i-4)...
                        - b_V*V_T...
                        - nu*V_T...
                        + nu/(N_A)*sum(V_R)...
                        + nu/(N_A)*sum(V_1)...
                        + nu/(N_A)*sum(V_2);
    else
        dydt(5*i-2) =   a_V*V_T*y(5*i-4)...
                        - b_V*V_T;
    end

    % VISITING HONEYBEES (regular/tolerant)
    % again, if there are no resident bees then none will visit
    if H_T > threshold
        dydt(5*i-1) =   - q_H*H_T*y(5*i-1)...
                        - alpha_H*Burden*y(5*i-1)...
                        - b_H*y(5*i-1)...
                        - nu*y(5*i-1)...
                        + nu/(N_A)*sum(H_1);
    else
        dydt(5*i-1) =   - q_H*H_T*y(5*i-1)...
                        - alpha_H*Burden*y(5*i-1)...
                        - b_H*y(5*i-1);
    end

    % VISITING HONEYBEES (tolerant/resistant)
    % again, if there are no resident bees then none will visit#
    if H_T > threshold
        dydt(5*i-0) =   - q_H*H_T*y(5*i-0)...
                        - alpha_H*Burden*y(5*i-0)...
                        - b_H*y(5*i-0)...
                        - nu*y(5*i-0)...
                        + nu/(N_A)*sum(H_2);
    else
        dydt(5*i-0) =   - q_H*H_T*y(5*i-0)...
                        - alpha_H*Burden*y(5*i-0)...
                        - b_H*y(5*i-0);
    end
end

end