clear all;
x0 = [2.5, 25, 7];
lb = [0.1 0.1 0.1];
start_u = -0.333;
start_y = -0.08839;
end_y = 1.5;
kk=400;

[nastawy_PID_fmincon,E]=fmincon((@(parameters) PID_fun(parameters,start_y,end_y,start_u)),x0,[],[],[],[],lb,[]);
%[nastawy_PID_ga, ~, ~] = ga((@(parameters) PID_fun(parameters)), 3, [], [], [], [], [0 0.1 0], []);
%save('optymalne_parametry_PID.mat', 'nastawy_PID_fmincon');

load('temp_to_del.mat');

% %%%ZAPIS DO PLIKU%%%%%%%
% odkomentuj by w³¹czyæ
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
%    disp('Nie udalo sie stworzyc folderów')
% end
% warning('on','all') %wlaczenie warningow

%wyniki symulacji
% figure; stairs(U);
% title('Sterowanie regulatora PID'); xlabel('k'); ylabel('wartoœæ sygna³u');
figure; stairs(Y); hold on; stairs(yzad,':'); legend('wyjœcie y(k)','wartoœæ zadana');
title('Wyjœcie regulatora PID'); xlabel('k'); ylabel('wartoœæ sygna³u');
