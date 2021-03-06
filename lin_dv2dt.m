function [outputArg1] = lin_dv2dt(V1,V2,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin)
     h1 = (V1/C1-h1lin^2)/(2*h1lin)+h1lin;
%     if h1 <0
%         h1 = 0;
%     end
     h2 = (V2/C2-h2lin^3)/(3*h2lin^2)+h2lin;
%     if h2 <0
%         h2 = 0;
%     end
  
    a = alfa1*(h1lin^(1/2)+1/(2*h1lin^(1/2))*(h1-h1lin));
%     if a < 0
%         a = 0;
%     end
    b = alfa2*(h2lin^(1/2)+1/(2*h2lin^(1/2))*(h2-h2lin));
%     if b < 0
%         b = 0;
%     end
    outputArg1 = a - b;
end

%     if alfa1/C1^(1/4)*(v1zero^(1/4)+(v1zero^(1/4)*(V1-v1zero)/(4*v1zero))) < 0
%         a = 0;
%     else
%         a = alfa1/C1^(1/4)*(v1zero^(1/4)+(v1zero^(1/4)*(V1-v1zero)/(4*v1zero)));
%     end
%     if alfa2/C2^(1/6)*(v2zero^(1/6)+(v2zero^(1/6)*(V2-v2zero)/(6*v2zero))) < 0
%         b = 0;
%     else
%         b = alfa2/C2^(1/6)*(v2zero^(1/6)+(v2zero^(1/6)*(V2-v2zero)/(6*v2zero))) < 0;
%     end
%     if a - b <0
%         outputArg1 = 0;
%     else
%         outputArg1 = a -b;
%     end
% % outputArg1 = alfa1*(V1/C1)^(1/4) - alfa2*(V2/C2)^(1/6);
% end