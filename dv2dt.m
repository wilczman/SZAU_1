function [outputArg1] = dv2dt(V1,V2,F1,Fd,alfa1,alfa2,C1,C2)
    a = alfa1*(V1/C1)^(1/4);
    if V1/C1 <0
        a = 0;
    end
%     if alfa1*(V1/C1)^(1/4) < 0
%         a = 0;
%     else
%         a = alfa1*(V1/C1)^(1/4);
%     end
    b = alfa2*(V2/C2)^(1/6);
    if V2/C2 <0
        b = 0;
    end
%     if alfa2*(V2/C2)^(1/6) < 0
%         b = 0;
%     else
%         b = alfa2*(V2/C2)^(1/6);
%     end
    outputArg1 = a - b;
%     if a - b <0
%         outputArg1 = 0;
%     else
%         outputArg1 = a -b;
%     end
% outputArg1 = alfa1*(V1/C1)^(1/4) - alfa2*(V2/C2)^(1/6);
end