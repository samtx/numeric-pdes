phi0{1} = @(s) -27/48*(s^3 - s^2 - 1/9*s + 1/9);
phi0{2} = @(s) 27/16*(s^3 - 1/3*s^2 - s + 1/3);
phi0{3} = @(s) -27/16*(s^3 + 1/3*s^2 - s - 1/3);
phi0{4} = @(s) 27/48*(s^3 + s^2 - 1/9*s - 1/9);

xx = linspace(-1,1)';
y = zeros(length(xx),4);
for i = 1:4
for j = 1:length(xx)
x = xx(j);
y(j,i) = phi0{i}(x);
end
end
x = [xx,xx,xx,xx];
plot(x,y)
