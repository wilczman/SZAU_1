%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                              %
%%    Przyk³adowy skrypt z implementacj¹ prostego regulatora DMC w wersji       %
%%  analitycznej i numerycznej.                                                 %
%%    Przyjête oznaczenia:                                                      %
%%    s - rzêdne odpowiedzi skokowej obiektu                                    %
%%    D - horyzont dynamiki                                                     %
%%    N - horyzont predykcji                                                    %
%%    Nu - horyzont sterowania                                                  %
%%    lambda - wspó³czynnik kary za przyrosty sterowania                        %
%%    Ke, Ku - parametry regulatora                                             %
%%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Piotr Marusak; 2001-2020
%   Wersja przeznaczona dla s³uchaczy wyk³adu Sztuczna inteligencja w automatyce
% (SZAU) prowadzonego na Wydziale Elektroniki i Technik Informacyjnych 
% Politechniki Warszawskiej.

clear all;

y_zad=1;
czas_sym=11;
fl_anal=1;

% Wektor z rzêdnymi odpowiedzi skokowej obiektu

s=[0.6000    0.8400    0.9360    0.9744    0.9898    0.9959    0.9984    0.9993    0.9997    0.9999];

% Horyzonty

D=6;
N=6;
Nu=3;

% Ograniczenie

umax=1.1;

% Wspó³czynnik kary za przyrosty sterowania

lambda=.1;

% Warunki pocz¹tkowe

y_pocz=0;
u_pocz=0;
yzad_pocz=0;

yk=y_pocz;
uk=u_pocz;
for i=1:D-1
   deltaupk(i)=0;
end

% Generowanie macierzy

M=zeros(N,Nu);
for i=1:N
   for j=1:Nu
      if (i>=j)
         M(i,j)=s(i-j+1);
      end;
   end;
end;

MP=zeros(N,D-1);
for i=1:N
   for j=1:D-1
      if i+j<=D
         MP(i,j)=s(i+j)-s(j);
      else
         MP(i,j)=s(D)-s(j);
      end;      
   end;
end;

% Obliczanie parametrów regulatora

I=eye(Nu);
K=((M'*M+lambda*I)^-1)*M';
Ku=K(1,:)*MP;
Ke=sum(K(1,:));

% Macierze dla QDMC

H=2*(M'*M+lambda*I);
A=zeros(Nu,Nu);
for i=1:Nu
   for j=1:Nu
      if (i>=j)
         A(i,j)=1;
      end
   end
end

% G³ówna pêtla programu
for i=1:czas_sym
   
   % Obiekt
   
   ukob=uk;
   if ukob>umax
      ukob=umax;
   end
   yk=0.4*yk+0.6*ukob;
   
   % Regulator
   
   ek=y_zad-yk;
   
   if (fl_anal==1)
      deltauk=Ke*ek-Ku*deltaupk';
      if uk+deltauk>umax
         deltauk=umax-uk;
      end   
   else
      f=-2*M'*(ek*ones(N,1)-MP*deltaupk');
      b=(umax-uk)*ones(Nu,1);
      x=quadprog(H,f,A,b);
      deltauk=x(1);
   end
   for n=D-1:-1:2;
      deltaupk(n)=deltaupk(n-1);
   end
   deltaupk(1)=deltauk;
   uk=uk+deltaupk(1);
   
   wyu(i)=uk;
   wyy(i)=yk;
end

% Graficzna prezentacja wyników obliczeñ
kolor='b';
figure;
stairs(0:czas_sym,[u_pocz wyu],kolor);hold on; grid on;
xlabel('czas [Tp]');ylabel('u')
axis([0 10 0 1.2]);
figure;
stairs(0:czas_sym,[yzad_pocz y_zad*ones(1,czas_sym)],'k--');
hold on; grid on;
plot(1:czas_sym,[wyy],kolor);
axis([0 10 0 1.2]);
xlabel('czas [Tp]');ylabel('y, y_{zad}')
