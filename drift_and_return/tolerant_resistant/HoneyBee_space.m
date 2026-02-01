function dydt = HoneyBee_space(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params)

N = length(y)/3;        % number of hives
% type - vector of length N, detailing what category of bee is in each hive
% if type(i) == 1, regular hive
% if type(i) == 2, tolerant hive
% if type(i) == 3, resistant hive

dydt = zeros(N*3, 1);

%% Set Parameters - fixed regardless of colony/hive type
% BB_p - Bee Brood Parameters
% HB_p - Honeybee Parameters
% VM_p - Varroa mite parameters
% SP_p - Space parameters

a_H = BB_p(1);      omega = BB_p(2);    b_B = BB_p(3);      m = BB_p(4);
q_H = HB_p(1);
b_V = VM_p(1);      alpha_B = VM_p(2);  sigma = VM_p(3);    
nu  = SP_p(1);      H_0 = SP_p(2);

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

H = y(2:3:end)';
V = y(3:3:end)';
Burden = V./H;
Burden(isinf(Burden)) = 0;

b_H = zeros(size(H));
b_H(H<10000) = HB_p(2);

for i = 1:N

    if types(i) == 1                     % Regular Bees
        a_V = type_params(1,1);
        alpha_H = type_params(1,2);
    elseif types(i) == 2                 % Tolerant Bees
        a_V = type_params(1,1);
        alpha_H = type_params(2,2);
    elseif types(i) == 3                 % Resistant Bees
        a_V = type_params(2,1);
        alpha_H = type_params(1,2);
    end

    % H(i) - current hive Honeybee density
    % V(i) - current hive mite density
    % Burden(i) - current burden

    visiting_dispersal  = nu/(N-1)*H(i)...
                            *sum(H./(H_0 + H)...
                            .*(Burden(i) - Burden));
    visitor_dispersal   = nu/(N-1)*H(i)/(H_0 + H(i))...
                            *sum(H ...
                            .*(Burden - Burden(i)));

    dydt(3*i-2) = a_H*y(3*i-1)^2/(omega^2 + y(3*i-1)^2) - b_B*y(3*i-2) - m*y(3*i-2) - alpha_B*y(3*i-0)*y(3*i-2)/(sigma + y(3*i-2));
    dydt(3*i-1) = m*y(3*i-2) - q_H*(y(3*i-1)^2) - alpha_H*y(3*i-0) - b_H(i)*y(3*i-1);
    dydt(3*i-0) = a_V*y(3*i-2)*y(3*i-0) - b_V*y(3*i-0) - visiting_dispersal + visitor_dispersal;
end

end