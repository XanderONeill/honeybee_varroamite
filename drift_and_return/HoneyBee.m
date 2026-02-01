function dydt = HoneyBee(t, y, BB_p, HB_p, VM_p)
dydt = zeros(3, 1);  
%   1 equation for brood, 1 for HB, 1 for varroa
% BB_p - all brood related parameters
% HB_p - all honeybeee related parameters
% VM_P - all varroa related parameters

a_H = BB_p(1);      omega = BB_p(2);    b_B = BB_p(3);      m = BB_p(4);
q_H = HB_p(1);      b_H = HB_p(2);
a_V = VM_p(1);      b_V = VM_p(2);      
alpha_B = VM_p(3);  sigma = VM_p(4);    alpha_H = VM_p(5);

% seasonal birth function, note t = 3-6 (summer) not required as this is
% set to a_H automatically. Cannot have final else, otherwise summer would
% have winter birth rate.
% if rem(t, 12) <= 3
%     a_H = 8*a_H/9;
% elseif (rem(t,12) > 6) && (rem(t,12) <= 9)
%     a_H = 6*a_H/9;
% elseif rem(t,12) > 9
%     a_H = 1*a_H/9;
% end

if rem(t, 12) <= 3
    a_H = 1*a_H/3;
elseif (rem(t,12) > 6) && (rem(t,12) <= 9)
    a_H = 1*a_H/3;
elseif rem(t,12) > 9
    a_H = 0*a_H/9;
end

if y(2) < 10000
    b_H = b_H; %#ok<ASGSL>
else
    b_H = 0;
end

dydt(1) = a_H*y(2)^2/(omega^2 + y(2)^2) - b_B*y(1) - m*y(1) - alpha_B*y(3)*y(1)/(sigma + y(1));
dydt(2) = m*y(1) - q_H*(y(2)^2) - alpha_H*y(3) - b_H*y(2);
dydt(3) = a_V*y(1)*y(3) - b_V*y(3);

end