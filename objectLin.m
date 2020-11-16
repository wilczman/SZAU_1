function [V1lin_out, V2lin_out] = objectLin(t, h,V1lin,V2lin,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin)
    
    k11=lin_dv1dt(V1lin(t),V2lin(t),F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin); %obliczanie wspó³czynników k dla obu zmiennych 
    k12=lin_dv2dt(V1lin(t),V2lin(t),F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k21=lin_dv1dt(V1lin(t)+0.5*h*k11,V2lin(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k22=lin_dv2dt(V1lin(t)+0.5*h*k11,V2lin(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k31=lin_dv1dt(V1lin(t)+0.5*h*k21,V2lin(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k32=lin_dv2dt(V1lin(t)+0.5*h*k21,V2lin(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k41=lin_dv1dt(V1lin(t)+h*k31,V2lin(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k42=lin_dv2dt(V1lin(t)+h*k31,V2lin(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    V1lin_out=V1lin(t)+1/6*h*(k11+2*k21+2*k31+k41);
    if V1lin_out<0
        V1lin_out = 0;
    end
    V2lin_out=V2lin(t)+1/6*h*(k12+2*k22+2*k32+k42);   
    if V2lin_out<0
        V2lin_out = 0;
    end   
end