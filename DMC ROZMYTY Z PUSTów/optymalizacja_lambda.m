clear all
liczba_regulatorow=5;
k=10;  %ile razy ma sie wykonac ga
nastawy_table=zeros(k,liczba_regulatorow);
error_table=zeros(k,1);

lb=zeros(1, liczba_regulatorow);
ub=zeros(1, liczba_regulatorow);
x0=zeros(1, liczba_regulatorow);
for i=1:liczba_regulatorow
    lb(i)=1;
    ub(i)=30;
    x0(i)=1;
end

save('optym_lambda_ws.mat', 'liczba_regulatorow')

options = gaoptimset('StallGenLimit', 10, 'PopulationSize', 30);
% for i=1:k
%     [nastawy_DMC_lambda, fval] = ga(@DMC_lambda_fun, liczba_regulatorow, [], [], [], [], lb, ub, [], int_vars, options);
%     
%     nastawy_table(i,:)=nastawy_DMC_lambda;   
%     error_table(i)=fval;
% end

[nastawy_DMC_lambda, fval] = fmincon((@(parameters) DMC_lambda_fun(parameters, liczba_regulatorow)), x0, [], [], [], [], lb, ub);

[min_val, indeks]=min(error_table);
DMC_lambda_fun(nastawy_DMC_lambda, liczba_regulatorow);
%nastawy_DMC_ga=round(nastawy_DMC_ga);
load('temp_lambda_DMC.mat');
figure; stairs(U);hold on
title('Sygna³ sterowania DMC'); xlabel('k');
figure; stairs(Y); hold on; 
stairs(yzad,':');     
title('Wyjœcie regulatora DMC'); xlabel('k'); ylabel('wartoœæ sygna³u');
legend('wyjœcie y(k)','wartoœæ zadana');


j=0; k=0;
for k=1:size(U,2)       
       j=j+1;
       table_U(1,k)=double(j);
       table_Y(1,k)=double(j);
       table_yzad(1,k)=double(j);
end
%macierze wykorzystane do zapisu
table_U(2,:)=U;
table_Y(2,:)=Y;
table_yzad(2,:)=yzad;

fname_U = sprintf('U.txt');
fname_Y = sprintf('Y.txt');
fname_yzad = sprintf('Yzad.txt');
fname_params = sprintf('params(%d).txt',liczba_regulatorow);
mkdir_status_U=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow));
mkdir_status_Y=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow));
mkdir_status_yzad=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow));
mkdir_status_params=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow));
if mkdir_status_U && mkdir_status_Y && mkdir_status_yzad && mkdir_status_params
    savdir_U =  sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow);
    savdir_Y = sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow);
    savdir_yzad = sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow);
    savdir_params = sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_lambda\\l_reg(%.0f)\\',liczba_regulatorow);
    
   fileID = fopen([savdir_U fname_U],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_U);
   fclose(fileID);
   fileID = fopen([savdir_Y fname_Y],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_Y);
   fclose(fileID);
   fileID = fopen([savdir_yzad fname_yzad],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_yzad);
   fileID = fopen([savdir_params fname_params],'w');
   fprintf(fileID,'lambdy= ');
   fprintf(fileID,'%f \t',nastawy_DMC_lambda);
   fprintf(fileID,'\r\n');
   fprintf(fileID,'E=%f ', fval);
   fclose(fileID);
else 
   disp('Nie udalo sie stworzyc folderów')
end
warning('on','all') %wlaczenie warningow
nastawy_DMC_lambda