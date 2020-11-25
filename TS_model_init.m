h2_pocz=15;
h2_koniec=45;
x = (h2_pocz:0.1:h2_koniec)';
kk=101;
liczba_regulatorow=5;

sigma=[3.2, 3.6, 2.9, 3.5, 3.5];
centra=linspace(h2_pocz,h2_koniec,liczba_regulatorow);
centra(3)=28.63;
centra(1)=centra(1)+2;
%centra(2)=centra(2);
save('TS_param.mat', 'centra','sigma')

