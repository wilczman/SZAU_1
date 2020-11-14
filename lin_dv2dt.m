function [outputArg1] = lin_dv2dt(V1,V2,F1,Fd,alfa1,alfa2,C1,C2,v1zero,v2zero)
    if alfa1/C1^(1/4)*(v1zero^(1/4)+(v1zero^(1/4)*(V1-v1zero)/(4*v1zero))) < 0
        a = 0;
    else
        a = alfa1/C1^(1/4)*(v1zero^(1/4)+(v1zero^(1/4)*(V1-v1zero)/(4*v1zero)));
    end
    if alfa2/C2^(1/6)*(v2zero^(1/6)+(v2zero^(1/6)*(V2-v2zero)/(6*v2zero))) < 0
        b = 0;
    else
        b = alfa2/C2^(1/6)*(v2zero^(1/6)+(v2zero^(1/6)*(V2-v2zero)/(6*v2zero))) < 0;
    end
    if a - b <0
        outputArg1 = 0;
    else
        outputArg1 = a -b;
    end
% outputArg1 = alfa1*(V1/C1)^(1/4) - alfa2*(V2/C2)^(1/6);
end