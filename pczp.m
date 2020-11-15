a(1:1818) = 50*sin(x/1818*2*pi)+100
b = repmat(a,1,11)
plot(b)
ylabel('x [µm]')
xlabel('Czas [µs]')