clear all;

%%%%%%%%%%% inicjalizacja %%%%%%%%%%
liczba_regulatorow = 5 ; %iloœæ rozmytych
 centra = linspace(-1,1,liczba_regulatorow);
%centra = [-1,1];
params_PID{1} = [ 0.734  2.398   1.668 ];
params_PID{2} = [ 9.784  3.220   0.722 ;0.100  2.239   2.580 ];
params_PID{3} = [9.784,3.220,0.722 ;0.734,2.398,1.668; 0.100  2.239   2.580 ];
params_PID{4} = [9.784,3.220,0.722 ; 5.238  2.218  1.532 ; 0.342  3.151  2.150; 0.100  2.239   2.580 ] ;
params_PID{5} = [9.784,3.220,0.722 ; 6.335  2.060 1.606 ; 0.734  2.398   1.668; 0.187  6.664  0.198 ;0.100  2.239   2.580 ];

%Punkt Pracy
Upp=0;
Ypp=0;

%Ograniczenia
u_min=-1;
u_max=1;
u_max=u_max-Upp;
u_min=u_min-Upp;
delta_u_max=0.1;

kk=1450; %d³ugoœæ symulacji

%deklaracja wektorów sygna³ów oraz b³êdów
for nr=1:liczba_regulatorow
    u{nr}=zeros(1,kk); 
end
U=zeros(1, kk);
U(:)=Upp;
Y=zeros(1, kk);
y=zeros(1, kk);
Y(1:12)=Ypp;
e=zeros(1, kk);
u_final=zeros(1,kk);
yzad=zeros(1, kk);

yzad(round(1*kk/7):round(2*kk/7))=4;
yzad(round(2*kk/7):round(3*kk/7))=0;
yzad(round(3*kk/7):round(4*kk/7))=0.2;
yzad(round(4*kk/7):round(5*kk/7))=1.1;
yzad(round(5*kk/7):round(7*kk/7))=3;

yzad=yzad-Ypp;

% %%%%%% TEST %%%%%%
% for nr = 1:liczba_regulatorow
%     x=linspace(-1,1,1000);
%     gs{nr}=gaussmf(x, [sigma(liczba_regulatorow),centra(nr)]);
%     hold on;
%     plot(gs{nr});
% end

load('optymalne_parametry_PID.mat');
% K=K_kryt*0.5; Ti=22; Td=7; %NAJLEPSZE EKSPERYMENTALNE PARAMETRY

for nr=1:liczba_regulatorow
    K(nr)=params_PID{liczba_regulatorow}(nr,1); 
    Ti(nr)=params_PID{liczba_regulatorow}(nr,2);
    Td(nr)=params_PID{liczba_regulatorow}(nr,3);  %parametry dobrane funkcj¹ fmincon
end

T=0.5;


for k=12:kk %g³ówna pêtla symulacyjna
    %symulacja obiektu    
    Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));
    y(k)=Y(k)-Ypp;
    %uchyb regulacji
    e(k)=yzad(k)-y(k);
    %sygna³ steruj¹cy regulatora PID
    %u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    
    for nr=1:liczba_regulatorow
        r2(nr)=K(nr)*Td(nr)/T;
        r1(nr)=K(nr)*(T/(2*Ti(nr))-2*Td(nr)/T-1);
        r0(nr)=K(nr)*(1+T/(2*Ti(nr))+Td(nr)/T);
        delta_u=r2(nr)*e(k-2)+r1(nr)*e(k-1)+r0(nr)*e(k);  
        u{nr}(k)=delta_u+ u_final(k-1);
    end
    
    
    %%%%%%------ROZMYCIE------%%%%%%%%
    for nr = 1:liczba_regulatorow
%         w(nr)=gs{nr}( round( (u_final(k-1)+1) /2 *999)+1);
        w(nr) = gaussmf(u_final(k-1), [gausy(liczba_regulatorow) centra(nr)]);
    end
    
    %œrednia wa¿ona z f.przynale¿noœci
    for nr=1:liczba_regulatorow 
        u_final(k) = u_final(k) + u{nr}(k) * w(nr); 
    end
    
    u_final(k) = u_final(k) / sum(w);
    
    %ograniczenie ró¿nic sygna³u steruj¹cego 
    if u_final(k)-u_final(k-1)>delta_u_max
         u_final(k)=u_final(k-1)+delta_u_max; 
    elseif u_final(k)-u_final(k-1)<(-delta_u_max)
         u_final(k)=u_final(k-1)-delta_u_max;
    end
    %ograniczenie wielkosci sygnalu
    if u_final(k)>u_max
        u_final(k)=u_max;
    elseif u_final(k)<u_min
        u_final(k)=u_min;
    end
    U(k)= u_final(k)+Upp;
     
end
yzad=yzad+Ypp;
wskaznik_jakosci=sum(e.^2);
%wyniki symulacji
figure; stairs(U);
title('Sterowanie regulatora PID'); xlabel('k'); ylabel('wartoœæ sygna³u');
figure; stairs(Y); hold on; stairs(yzad,':'); legend('wyjœcie y(k)','wartoœæ zadana');
title('Wyjœcie regulatora PID'); xlabel('k'); ylabel('wartoœæ sygna³u');




%%%ZAPIS DO PLIKU%%%%%%%
% % %   odkomentuj by w³¹czyæ
j=0; k=0;
for k=1:kk       
       j=j+1;
       table_U(1,k)=double(j);
       table_Y(1,k)=double(j);
       table_yzad(1,k)=double(j);
end
%macierze wykorzystane do zapisu
table_U(2,:)=U;
table_Y(2,:)=Y;
table_yzad(2,:)=yzad;

fname_U = sprintf('outU_liczba_regulatorow(%d).txt',liczba_regulatorow);
fname_Y = sprintf('outY_liczba_regulatorow(%d).txt',liczba_regulatorow);
fname_yzad = sprintf('outYzad_liczba_regulatorow(%d).txt',liczba_regulatorow);
mkdir_status_U=mkdir(sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\wyniki_PID_rozmyty'));
if mkdir_status_U 
   savdir_U =  sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\wyniki_PID_rozmyty\\');
    
   fileID = fopen([savdir_U fname_U],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_U);
   fclose(fileID);
   fileID = fopen([savdir_U fname_Y],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_Y);
   fclose(fileID);
   fileID = fopen([savdir_U fname_yzad],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_yzad);
   fclose(fileID);
else 
   disp('Nie udalo sie stworzyc folderów')
end
warning('on','all') %wlaczenie warningow