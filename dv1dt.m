function [outputArg1] = dv1dt(V1,V2,F1,Fd,alfa1,alfa2,C1,C2)
outputArg1 = F1 + Fd - alfa1*(V1/C1)^(1/4);
end