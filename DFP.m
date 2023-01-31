m = 5;
n = 3;
y_0 = 1;
w = sym('w',[1 n]);
v = sym('v',[1 n]);
theta = sym('theta',[1 n]);
s1 = 0;
s2 = 0;
y = 0;
for i = 1:m
    for j = 1:n
        s1 = s1 + v(j) * 1 / (1 + exp(theta(j) - w(j) * i / m));
        s2 = s1+v(j)*w(j)*exp(theta(j)-w(j)*i/m)/(1+exp(theta(j)-w(j)*i/m))^2;
    end
    y = y+((y_0+i/m*s1-(i/m)^2)*(i/m+(1+3*(i/m)^2)/(1+i/m+(i/m)^3))-2*i/m+s1+i/m*s2)^2;
end
f(w,v,theta) = y;
p = [w,v,theta];
g = sym('g',[3*n 1]);
for i = 1:3*n
    g(i) = diff(f,p(i));
end
p1 = ones(3*n,1);
H = eye(3*n);
epsilon = 1.5;
k = 0;
while(f(p1(1),p1(2),p1(3),p1(4),p1(5),p1(6),p1(7),p1(8),p1(9))>epsilon && k < 1000)
    g1 = zeros(3*n,1);
    for i = 1:3*n
        G(p) = g(i);
        g1(i) = vpa(G(p1(1),p1(2),p1(3),p1(4),p1(5),p1(6),p1(7),p1(8),p1(9)));
    end
    d = - H * g1;
    p2 = p1 + 0.0001 * d;
    g2 = zeros(3*n,1);
    for i = 1:3*n
        G(p) = g(i);
        g2 = vpa(G(p2(1),p2(2),p2(3),p2(4),p2(5),p2(6),p1(7),p1(8),p1(9)));
    end
    s = d;
    y = g2 - g1;
    T1 = s * s.';
    T2 = s.' * y;
    D1 = T1 / T2;
    T3 = H * y;
    T4 = T3 * T3.';
    T5 = y.' * T3;
    D2 = T4 / T5;
    H = H + D1 - D2;
    p1 = p2
    k = k + 1
    vpa(f(p1(1),p1(2),p1(3),p1(4),p1(5),p1(6),p1(7),p1(8),p1(9)))
end