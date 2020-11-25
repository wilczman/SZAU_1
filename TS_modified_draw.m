clear all;

liczba_regulatorow = 5;
%zakres wysokosci
figure;
hold on;

h2_pocz=15;
h2_koniec=45;
x = (h2_pocz:0.1:h2_koniec)';
kk=101;
center=linspace(h2_pocz,h2_koniec,liczba_regulatorow);
load('TS_param.mat')


for nr=1:liczba_regulatorow
    y{nr}=gaussmf(x, [sigma(nr) centra(nr)]);
    
    %subplot(4,1,i);
    str=sprintf('Liczba modeli = %d', liczba_regulatorow);
    title(str);
    xlabel('h_2')
    hold on
    plot(x, y{nr});
    
end
