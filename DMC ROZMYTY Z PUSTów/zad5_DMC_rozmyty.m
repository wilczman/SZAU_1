clear all

%%%%%%%    inicjalizacja    %%%%%%

liczba_regulatorow = 2; %ilo�� rozmytych
centra = linspace(-1,1,liczba_regulatorow);
% centra = [1];
params_DMC{1} = [ 29  10  5];
params_DMC{2} = [ 20 26 5;  85 12 6];
params_DMC{3} = [ 20 26 5; 29  10  5; 85 12 6 ];
params_DMC{4} = [ 20 26 5; 40 22 5 ; 60 8 71; 85 12 6 ] ;
params_DMC{5} = [ 20 26 5; 40 22 5 ; 29  10  5; 60 8 71; 85 12 6 ];

%%%%%%%%Punkt Pracy%%%%%%%%
Upp=0;
Ypp=0;

%%%%%%%%%Ograniczenia%%%%%%%
u_min=-1-Upp;
u_max=1-Upp;


%%%%%%%%%-----Parametry DMC------%%%%%%%
lambda=[1,1,1,1,1]; %parametr lambda np. 1
D=[50,50,50,50,50]; %horyzont dynamiki (D)
N=[50,50,50,50,50];%horyzont predykcji (N)
Nu=[30,30,30,30,30]; %horyzont sterowania (Nu)(ilosc przyszlych przyrostow wartosci sterowania)

<<<<<<< HEAD
kk=450+max(N);
=======
>>>>>>> b50d5f5d2b0370571bfdf8d1fbd682ad81a06cbf
%load('optymalne_parametry_DMC.mat');
%parametry dobrane przez funkcj� fmincon()
%D=150; N=nastawy_DMC_fmincon(1); Nu=nastawy_DMC_fmincon(2); lambda=nastawy_DMC_fmincon(3);
%D=150; N=35; Nu=5; lambda=1 %NAJLEPSZE EKSPERYMENTALNE PARAMETRY SWIATA
%D=200; N=50; Nu=10; % robocze parametry


%%%%%%deklaracja wektor�w sygna��w oraz b��d�w%%%%%%
for nr=1:liczba_regulatorow
    
    D(nr)=params_DMC{liczba_regulatorow}(nr,1); 
    N(nr)=params_DMC{liczba_regulatorow}(nr,2);
    Nu(nr)=params_DMC{liczba_regulatorow}(nr,3);  %parametry dobrane funkcj� ga

end
kk=450;


y=zeros(1, kk);
Y=zeros(1, kk); 
Y(1:12)=Ypp;
U=zeros(1, kk);  
U(:)=Upp;
e=zeros(1, kk);
u_final=zeros(1,kk);

<<<<<<< HEAD
yzad(round(1*kk/7):round(2*kk/7))=4;
yzad(round(2*kk/7):round(3*kk/7))=0;
yzad(round(3*kk/7):round(4*kk/7))=0.2;
yzad(round(4*kk/7):round(5*kk/7))=1.1;
yzad(round(5*kk/7):round(7*kk/7))=3;
=======
yzad=zeros(1, kk);
>>>>>>> b50d5f5d2b0370571bfdf8d1fbd682ad81a06cbf

yzad(round(1*kk/7):round(2*kk/7))=4;
yzad(round(2*kk/7):round(3*kk/7))=0;
yzad(round(3*kk/7):round(4*kk/7))=0.2;
yzad(round(4*kk/7):round(5*kk/7))=1.1;
yzad(round(5*kk/7):round(7*kk/7))=3;
yzad=yzad-Ypp;

%%%%%%%   tycie odpowiedzi skokowe  %%%%%%% 
for nr=1:liczba_regulatorow
    u{nr}=zeros(1,kk);
    %odp_skok{nr} = gotowa_odp_skokowa;
    odp_skok{nr} = fun_odp_skok(centra(nr),0.1);
end

%%%%%%% Macierze M, K, Mp %%%%%%%
for nr=1:liczba_regulatorow
    Mp{nr}=zeros(N(nr),D(nr)-1);        %macierz ma wymiary Nx(D-1)
    for i=1:D(nr)-1 %wypelnianie macierzy Mp
       Mp{nr}(1:N(nr), i)=odp_skok{nr}(i+1:N(nr)+i)-odp_skok{nr}(i);
    end
    %macierz wsp�czynnik�w odpowiedzi skokowej wymiary(NxNu)
    M=zeros(N(nr), Nu(nr));  
    i=0;
    for j=1:Nu(nr)  %wypelnianie macierzy trojkatnej dolnej M
       M(j:N(nr),j)=odp_skok{nr}(1:N(nr)-i).';  
       i=i+1;
    end

    I=eye(Nu(nr));  %tworzenie macierzy jednostkowej o wymiarach NuxNu
    K{nr}=inv(M.'*M+lambda(nr)*I)*M.';   %macierz K

    deltaUP{nr}(1:D(nr)-1,1)=0;
    deltaU{nr}=0;
end

%%%%%%%%% Algorytm DMC %%%%%%%%%
for k=12:kk-max(N) 
    %symulacja obiektu i regulatora
    Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));
    y(k)=Y(k)-Ypp;
    e(k)=yzad(k)-y(k);
    for nr=1:liczba_regulatorow
        deltaUP{nr}(2:D(nr)-1)=deltaUP{nr}(1:D(nr)-2);
        deltaUP{nr}(1) = u_final(k-1)-u_final(k-2);
        Y0{nr}=Mp{nr}*deltaUP{nr}+y(k);
        Yzad=yzad(k+1:k+N(nr));
        deltaU{nr}=K{nr}*(Yzad-Y0{nr});	
        delta_u{nr}=deltaU{nr}(1);

        u{nr}(k)=u_final(k-1)+delta_u{nr};
    end
    
    %%%%%%%------ROZMYCIE------%%%%%%%%
    for nr = 1:liczba_regulatorow
        w(nr) = gaussmf(u_final(k-1), [gausy(liczba_regulatorow) centra(nr)]);
    end
    
    %�rednia wa�ona z f.przynale�no�ci
    for nr=1:liczba_regulatorow 
        u_final(k) = u_final(k) + u{nr}(k) * w(nr); 
    end
    u_final(k) = u_final(k) / sum(w);
    
    %ograniczenie sygna�u steruj�cego
    if u_final(k)>u_max
        u_final(k)=u_max;
    elseif u_final(k)<u_min
        u_final(k)=u_min;
    end
     U(k)=u_final(k)+Upp;
     
end
wskaznik_jakosci=sum(e.^2);
yzad=yzad(1:kk-N(nr))+Ypp;
Y=Y(1:kk-N(nr));
U=U(1:kk-N(nr));

% for nr=1:liczba_regulatorow
%     wskaznik_jakosci=sum(e.^2);
%     yzad=yzad(1:kk-N(nr))+Ypp;
%     Y=Y(1:kk-N(nr));
%     U=U(1:kk-N(nr));
% end

    %%%%%%%%prezentacja wynik�w symulacji%%%%%%%%
figure;
stairs(U);hold on; xlim([0, kk-N(nr)]);
title('Sygna� sterowania DMC'); xlabel('k');ylabel('warto�� sygna�u');

figure;
stairs(Y); hold on; 
<<<<<<< HEAD
stairs(yzad,':'); xlim([0, kk]); %ylim([-1, 7]);
=======
stairs(yzad,':'); xlim([0, kk-N(nr)]); 
>>>>>>> b50d5f5d2b0370571bfdf8d1fbd682ad81a06cbf
title('Wyj�cie regulatora DMC'); xlabel('k'); ylabel('warto�� sygna�u');
legend('wyj�cie y(k)','warto�� zadana');

% %%%ZAPIS DO PLIKU%%%%%%%
% odkomentuj by w��czy�
% j=0; k=0;
% for k=1:kk       
%        j=j+1;
%        table_U(1,k)=double(j);
%        table_Y(1,k)=double(j);
%        table_yzad(1,k)=double(j);
% end
% %macierze wykorzystane do zapisu
% table_U(2,:)=U;
% table_Y(2,:)=Y;
% table_yzad(2,:)=yzad;
% 
% fname_U = sprintf('outU_punkt(%.2f).txt',start_u);
% fname_Y = sprintf('outY_punkt(%.2f).txt',start_u);
% fname_yzad = sprintf('outYzad_punkt(%.2f).txt',start_u);
% fname_params = sprintf('params_punkt(%.2f).txt',start_u);
% mkdir_status_U=mkdir(sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\U\\',start_u));
% mkdir_status_Y=mkdir(sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\Y\\',start_u));
% mkdir_status_yzad=mkdir(sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\Yzad\\',start_u));
% mkdir_status_params=mkdir(sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\params\\',start_u));
% if mkdir_status_U && mkdir_status_Y && mkdir_status_yzad && mkdir_status_params
%     savdir_U =  sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\U\\',start_u);
%     savdir_Y = sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\Y\\',start_u);
%     savdir_yzad = sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\Yzad\\',start_u);
%     savdir_params = sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID_rozmyty\\fmincon\\punkt(%.2f)\\params\\',start_u);
%     
%    fileID = fopen([savdir_U fname_U],'w');
%    fprintf(fileID,'%6.3f %6.3f\r\n',table_U);
%    fclose(fileID);
%    fileID = fopen([savdir_Y fname_Y],'w');
%    fprintf(fileID,'%6.3f %6.3f\r\n',table_Y);
%    fclose(fileID);
%    fileID = fopen([savdir_yzad fname_yzad],'w');
%    fprintf(fileID,'%6.3f %6.3f\r\n',table_yzad);
%    fileID = fopen([savdir_params fname_params],'w');
%    fprintf(fileID,'%6.3f %6.3f %6.3f\r\n',nastawy_PID_fmincon);
%    fclose(fileID);
% else 
%    disp('Nie udalo sie stworzyc folder�w')
% end
% warning('on','all') %wlaczenie warningow
