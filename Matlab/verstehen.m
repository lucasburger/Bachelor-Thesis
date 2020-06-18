% verstehen der funktion aput aus günther jüngel, finanzderivate mit matlab

K = 100;
T = 1;
sigma = 0.4;
omega = 1.5;
r = 0.1 ;

a = 5; N = 10; M = 5;
T0 = sigma^2*T/2;           % transformierte Endzeit
h = 2*a/N;                  % Ortsschrittweite
s = T0/M;                   % Zeitschrittweite

% Abkürzungen
alpha = s/h^2; k = 2*r/sigma^2; x = [-a:h:a]; S = K*exp(x);
% Definition von f(i,j)
S
pause 

for j = 1:M+1
    f(:,j) = exp(0.5*(k-1)*x + 0.25*(k+1)^2*(j-1)*s*ones(1,N+1)).*max(0,1-exp(x));
end
u = f(:,1)'; loop = zeros(1,N+1);  % setzt M < N+2 voraus
f

pause
% Projektions-SOR-Verfahren
for j = 1:M
    u(1) = f(1,j+1); u(N+1) = f(N+1,j+1);
    b(2:N) = (alpha/2)*(u(1:N-1)+u(3:N+1)) + (1-alpha)*u(2:N);
    b(1) = (1-alpha)*u(1) + (alpha/2)*u(2) + (alpha/2)*(f(1,j)+f(1,j+1));
    b(N+1) = (alpha/2)*u(N) + (1-alpha)*u(N+1);
    tol = 1e-10; error = 1;
    while error > tol
        error = 0;
        for i = 2:N
             z = ((alpha/2)*(u(i-1)+u(i+1))+b(i))/(alpha+1);
             z = max(u(i) + omega*(z-u(i)),f(i,j+1)');
             error = error + (u(i)-z)^2;
             u(i) = z;
        end
        loop(j) = loop(j) + 1;
    end
    u
    u'==f(:,j)
    disp('fertig')
    pause
end

% Rücktransformation
P = K*exp(-0.25*(k+1)^2*T0)*exp(-0.5*(k-1)*x).*u; 
result = [S;P;loop];