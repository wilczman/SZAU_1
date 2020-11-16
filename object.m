function [V1_out, V2_out] = object(t, h ,V1,V2,F1,Fd,alfa1,alfa2,C1,C2)
    
    k11=dv1dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2); %obliczanie wspó³czynników k dla obu zmiennych 
    k12=dv2dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2);
    k21=dv1dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k22=dv2dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k31=dv1dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k32=dv2dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k41=dv1dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    k42=dv2dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    
    
    V1_out=V1(t)+1/6*h*(k11+2*k21+2*k31+k41);
    if V1_out<0
        V1_out = 0;
    end
    V2_out=V2(t)+1/6*h*(k12+2*k22+2*k32+k42);   
    if V2_out<0
        V2_out = 0;
    end
end