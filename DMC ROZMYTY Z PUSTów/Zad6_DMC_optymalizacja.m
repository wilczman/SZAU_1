clear all
start_u=1;
start_y=0;
end_y=1;
skok_odp=-0.1;
k=10;  %ile razy ma sie wykonac ga
nastawy_table=zeros(k,3);
error_table=zeros(k,1);


lb = [20 5 5];
ub = [90 90 90];

save('skoki.mat', 'start_u', 'start_y', 'end_y', 'skok_odp')

options = gaoptimset('StallGenLimit', 10, 'PopulationSize', 30);
for i=1:k
    [nastawy_DMC_ga, fval] = ga(@DMC_fun, 3, [], [], [], [], lb, ub, [], [1 2 3], options);
    nastawy_table(i,:)=nastawy_DMC_ga;   
    error_table(i)=fval;
end
[min_val, indeks]=min(error_table);
DMC_fun(nastawy_table(indeks, :));
%nastawy_DMC_ga=round(nastawy_DMC_ga);
load('temp_DMC.mat');
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
fname_params = sprintf('params_punkt(%.2f).txt',start_u);
mkdir_status_U=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\',start_u));
mkdir_status_Y=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\',start_u));
mkdir_status_yzad=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\',start_u));
mkdir_status_params=mkdir(sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\params\\',start_u));
if mkdir_status_U && mkdir_status_Y && mkdir_status_yzad && mkdir_status_params
    savdir_U =  sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\',start_u);
    savdir_Y = sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\',start_u);
    savdir_yzad = sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\',start_u);
    savdir_params = sprintf('D:\\OneDrive\\Documents\\GitHub\\PUST_3\\strojenie_DMC_rozmyty\\punkt(%.2f)\\params\\',start_u);
    
   fileID = fopen([savdir_U fname_U],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_U);
   fclose(fileID);
   fileID = fopen([savdir_Y fname_Y],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_Y);
   fclose(fileID);
   fileID = fopen([savdir_yzad fname_yzad],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_yzad);
   fileID = fopen([savdir_params fname_params],'w');
   fprintf(fileID,'D=%6.0f N=%6.0f Nu=%6.0f E=%f\r\n',nastawy_table(indeks,1),nastawy_table(indeks,2), nastawy_table(indeks,3), fval);
   fprintf(fileID,'start_u=%6.3f start_y=%6.3f end_y=%6.3f skok_odp=%f\r\n',start_u,start_y, end_y, skok_odp);
   fclose(fileID);
else 
   disp('Nie udalo sie stworzyc folderów')
end
warning('on','all') %wlaczenie warningow
nastawy_table(indeks,:)