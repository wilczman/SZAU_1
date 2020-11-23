clear all;
%% Sta³e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
h1zero = 24;
h2zero = 24;
h1lin = 28.63;
h2lin = 28.63;
timespan = 10000;
tau = 50;

%% I
Fd = 14.2; %dop³yw zak³ócaj¹cy
F1ster(1:4000) = 98.5; %dop³yw wody do zbiornika
h=1;  %ustawienie kroku
t=tau+1;
% starttime=1; %pocz¹tek przedzia³u%
warpoczv1=C1*h1lin^2; %warunki pocz¹tkowe
warpoczv2=C2*h2lin^3;
V1lin(1:timespan)=warpoczv1;
V2lin(1:timespan)=warpoczv2;

%% DMC
% load('odp_skokH28.mat')
load('rw_odp_skok.mat','s')
final_odp_skok=s;
D=5000; N=200; Nu=50;lambda = 1;delta_u_max = 200;u_max = 200; u_min = 0;
Upp=98.5;
Ypp=0;%24;
kk = timespan;
%%%%%%deklaracja wektorów sygna³ów oraz b³êdów%%%%%%
U=zeros(1, kk);
u=zeros(kk);
% U(:)=Upp;
Y = zeros(1,kk) + h2lin;
y = zeros(1,kk)+ h2lin;
% Y(1:tau)=Ypp;
yzad=zeros(1, kk);
yzad(1:kk)= 29;
e = zeros(1, kk);%yzad-Y;
% yzad(round(kk/5):round(2*kk/5))=2.1;
% yzad(round(2*kk/5):round(3*kk/5))=2;
% yzad(round(3*kk/5):round(4*kk/5))=2.3;
% yzad(round(4*kk/5):round(5*kk/5))=2.4;
% yzad=yzad-Ypp;

yzad(1:4000) = h2lin;
yzad(4000:7000) = h2lin+5;
yzad(7000:10000) = h2lin-3;

Mp=zeros(N,D-1);        %macierz ma wymiary Nx(D-1)
for i=1:D-1 %wypelnianie macierzy Mp
   Mp(1:N, i)=final_odp_skok(i+1:N+i)-final_odp_skok(i);
end
%macierz wspó³czynników odpowiedzi skokowej wymiary(NxNu)
M=zeros(N, Nu);  
i=0;
for j=1:Nu  %wypelnianie macierzy trojkatnej dolnej M
   M(j:N,j)=final_odp_skok(1:N-i).';  
   i=i+1;
end

I=eye(Nu);              %tworzenie macierzy jednostkowej o wymiarach NuxNu
K=inv(M.'*M+lambda*I)*M.';   %macierz K

deltaUP(1:D-1,1)=0;
deltaU=0;

%%%%%%%%% Algorytm DMC %%%%%%%%%
for t=tau+2:kk-N %symulacja obiektu i regulatora
    
    %###########
    %F1 = F1ster(t-tau);
    [V1lin(t), V2lin(t)] = objectLin(t-1,h,V1lin,V2lin,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    %############
    Y(t) = (V2lin(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin;%nthroot(V2lin(t)/C2, 3);
    
    y(t)=Y(t)-Ypp;
    e(t)=yzad(t)-y(t);
    deltaUP(2:D-1)=deltaUP(1:D-2);
    deltaUP(1) = u(t-1)-u(t-2);  
    Y0=Mp*deltaUP+y(t);
    Yzad=yzad(t)*ones(N,1);
    %Yzad=yzad(t+1:t+N);
    deltaU=K*(Yzad-Y0);	
    delta_u=deltaU(1);

    %ograniczenie ró¿nic sygna³u steruj¹cego 
%     if delta_u>delta_u_max
%          delta_u=delta_u_max; 
%     elseif delta_u<(-delta_u_max)
%          delta_u=-delta_u_max;
%     end

    u(t)=u(t-1)+delta_u;
    %ograniczenie sygna³u steruj¹cego
%     if u(t)>u_max
%         u(t)=u_max;
%     elseif u(t)<u_min
%         u(t)=u_min;
%     end
     U(t)=u(t)+Upp;
     F1ster(t) = U(t);   
end

figure;title('wykres F1ster');hold on
plot(F1ster(1:timespan-N));hold off
figure;title('wykres Y');hold on;
plot(Y(1:timespan-N)-Ypp);plot(yzad(1:timespan-N));legend('Y','yzad');ylim([0,40]);hold off;

% %% II
% F1in(1:timespan) = 0; %dop³yw wody do zbiornika - wielkoœæ steruj¹ca
% F1ster(1:2000) = 0; %dop³yw wody do zbiornika 
% F1ster(2000:timespan) = 98.5;
% Fd = 14.2;%14.2; %dop³yw zak³ócaj¹cy
% h=1;  %ustawienie kroku
% t=tau+1;
% % starttime=1; %pocz¹tek przedzia³u%
% warpoczv1=C1*h1zero^2; %warunki pocz¹tkowe
% warpoczv2=C2*h2zero^3;
% V1lin(1:timespan)=warpoczv1;
% V2lin(1:timespan)=warpoczv2;
% 
% %% symulacja
% while t<timespan %wykonuj na przedziale [0,15]
%     F1 = F1ster(t-tau);
%     k11=lin_dv1dt(V1lin(t),V2lin(t),F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin); %obliczanie wspó³czynników k dla obu zmiennych 
%     k12=lin_dv2dt(V1lin(t),V2lin(t),F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
%     k21=lin_dv1dt(V1lin(t)+0.5*h*k11,V2lin(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
%     k22=lin_dv2dt(V1lin(t)+0.5*h*k11,V2lin(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
%     k31=lin_dv1dt(V1lin(t)+0.5*h*k21,V2lin(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
%     k32=lin_dv2dt(V1lin(t)+0.5*h*k21,V2lin(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
%     k41=lin_dv1dt(V1lin(t)+h*k31,V2lin(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
%     k42=lin_dv2dt(V1lin(t)+h*k31,V2lin(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
%     V1lin(t+1)=V1lin(t)+1/6*h*(k11+2*k21+2*k31+k41);
%     if V1lin(t+1)<0
%         V1lin(t+1) = 0;
%     end
%     V2lin(t+1)=V2lin(t)+1/6*h*(k12+2*k22+2*k32+k42);   
%     if V2lin(t+1)<0
%         V2lin(t+1) = 0;
%     end
%     t=t+h;
% end
% 
% 
% %% rysowanie wykresu
% figure;hold on
% title('wykres V1lin,V2lin')
% plot(V1lin);
% plot(V2lin);
% legend('V1lin','V2lin');
% hold off;
% 
% figure;hold on
% title('wykres h1,h2')
% plot(nthroot(V1lin/C1, 2));
% plot(nthroot(V2lin/C2, 3));
% legend('h1','h2');
% hold off;


% figure
% odp_skok = nthroot(V2lin/C2, 3);
% odp_skok = odp_skok(2000:timespan-1)
% final_odp_skok = (odp_skok-28.77)/10
% plot(final_odp_skok)
% save odp_skokH28.mat