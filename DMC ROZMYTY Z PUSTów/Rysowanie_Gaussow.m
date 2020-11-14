clear all;

liczba_regulatorow = 5;
x = (-1:0.02:1)';
kk=101;
center=linspace(-1,1,liczba_regulatorow);

figure;
hold on;
for nr=1:liczba_regulatorow
    y{nr}=gaussmf(x, [gausy(liczba_regulatorow) center(nr)]);
    plot(x, y{nr});
end
        
for nr=1:liczba_regulatorow
    
           %%%ZAPIS DO PLIKU%%%%%%%
    j=0; k=0;
    for k=1:kk       
           j=j+1;
           table_Y(1,k)=double(j);
    end
    %macierze wykorzystane do zapisu
    table_Y(2,:)=y{nr};
    table_Y(1,:)=x;

    fname_Y = sprintf('Y_liczba_regulatorow_(%d)_nr_wykresu_(%d).txt',liczba_regulatorow,nr);
    mkdir_status_Y=mkdir(sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\wykres_gauss\\liczba_regulatorow_(%d)',liczba_regulatorow));
    if mkdir_status_Y
       savdir_Y = sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\wykres_gauss\\liczba_regulatorow_(%d)\\',liczba_regulatorow);

       fileID = fopen([savdir_Y fname_Y],'w');
       fprintf(fileID,'%6.3f %6.3f\r\n',table_Y);
       fclose(fileID);
    else 
       disp('Nie udalo sie stworzyc folderów')
    end
    warning('on','all') %wlaczenie warningow
end