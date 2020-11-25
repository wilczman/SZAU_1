function [c, ceq]=u_limit_fun(x)
    load('lastU.mat')
    c=-(x(1,1)+last_U);
    ceq=[];
end