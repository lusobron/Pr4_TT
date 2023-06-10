% Algoritmos de recursión necesarios para la rutina lambertbattin
function seebat = seebatt( v )

% Implementación
c(1) =    0.2;
c(2) =    9.0 /  35.0;
c(3) =   16.0 /  63.0;
c(4) =   25.0 /  99.0;
% ...
c(21)=  484.0 /1935.0;
sqrtopv= sqrt(1.0 + v);
eta    = v / ( 1.0 + sqrtopv )^2;

% Proceso hacia adelante
delold = 1.0;
termold= c(1);   % * eta
sum1   = termold;
i= 1; 
while ((i <= 20) && (abs(termold) > 0.00000001 ))
    del  = 1.0 / ( 1.0 + c(i+1)*eta*delold );
    term = termold * (del - 1.0);
    sum1 = sum1 + term;
    i    = i + 1;
    delold = del;
    termold= term;
end

% Cálculo final
seebat = 1.0/ ((1.0/(8.0*(1.0+sqrtopv))) * ( 3.0 + sum1 / ( 1.0+eta*sum1 ) ) );