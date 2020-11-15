function [outputArg1] = lin_dv1dt(V1,V2,F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero)
    h1 = (V1/C1-h1zero^2)/(2*h1zero)+h1zero;
    if h1 <0
        h1 = 0;
    end
  
    a =  F1 + Fd;
    if F1 + Fd < 0
        a = 0;
    end
    b = alfa1*(h1zero^(1/2)+1/(2*h1zero^(1/2))*(h1-h1zero));
    if b < 0
        b = 0;
    end
    outputArg1 = a -b;
    if a - b <0
        outputArg1 = 0;
    end
end

% function [outputArg1] = lin_dv1dt(V1,V2,F1,Fd,alfa1,alfa2,C1,C2,v1zero,v2zero)
%     if F1 + Fd < 0
%         a = 0;
%     else
%         a =  F1 + Fd;
%     end
%     if alfa1/C1^(1/4)*(v1zero^(1/4)+(v1zero^(1/4)*(V1-v1zero)/(4*v1zero))) < 0
%         b = 0;
%     else
%         b = alfa1/C1^(1/4)*(v1zero^(1/4)+(v1zero^(1/4)*(V1-v1zero)/(4*v1zero)));
%     end
%     if a - b <0
%         outputArg1 = 0;
%     else
%         outputArg1 = a -b;
%     end
% % outputArg1 = F1 + Fd - alfa1*(V1/C1)^(1/4);
% end