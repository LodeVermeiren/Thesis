%-----------------------------------------------------------------------
%           General constants needed for all approximations
%-----------------------------------------------------------------------

%Defined Catalan's constant using G = vpa(catalan,1000) which is precise enough
G = 0.9159655941772190150546035149323841107741493742816721342664981196217630197762547694793565129261151062485744226191961995790358988033258590594315947374811584069953320287733194605190387274781640878659090247064841521630002287276409423882599577415088163974702524820115607076448838078733704899008647751132259971343407485407553230768565335768095835260219382323950800720680355761048235733942319149829836189977069036404180862179411019175327431499782339761055122477953032487537187866582808236057022559419481809753509711315712615804242723636439850017382875977976530683700929808738874956108936597719409687268444416680462162433986483891628044828150627302274207388431172218272190472255870531908685735423498539498309919115967388464508615152499624237043745177737235177544070853846440132174839299994757244619975496197587064007474870701490937678873045869979860644874974643872062385137123927363049985035392239287879790633644032354784535851927777787270906083031994301332316712476158709792455479119092126201854803963934243;
%Defining the fucntion evaluation f_2(-1)
f2eval =-4*G;
%Setting the scope of degree n
limit = 15;

%-----------------------------------------------------------------------
%                   Type-II Hermite-Padé approximations
%-----------------------------------------------------------------------

%P will be the vector containing the polynomial evaluations P_n,n(-1)
P = zeros(1,limit);
%Pint will be the vector containing the polynomial evaluations P_n,n(-1)
%multiplied by the constant to make them integers
Pint = zeros(1,limit);
%R will be the vector containing the polynomial evaluations R_n,n(-1)
R = zeros(1,limit);
%Rint will be the vector containing the polynomial evaluations R_n,n(-1)
%multiplied by the constant to make them integers
Rint = zeros(1,limit);
for n = 1:limit
    nodds = sym(1:2:4*n-1);
    [P(n),R(n)] = typeII(n,n);
    Pint(n) = lcm(nodds)^2*4^(2*n)*P(n);
    Rint(n) = lcm(nodds)^2*4^(2*n)*R(n);
end

%convergence should go to zero; according to our results
convergence = f2eval*P-R;
%f2approxs are the approximations for -4G
f2approxs = R./P;
%approxs are the approximations for G
approxs = f2approxs/(-4);
%error is the difference between approxs and G
error = abs(approxs-G);

%We plot error logarithmically
figure
semilogy(1:limit,error)
title('Logarithmic plot of the Type-II error')
xlabel('n')
ylabel('error')

%We plot approxs against G
figure
plot(approxs)
hold on
plot(repmat(G,[1,limit]))
title("Type-II approximations vs Catalan's constant")
xlabel('n')
hold off

%intconvergence should NOT go to zero; according to our results
intconvergence = f2eval*Pint-Rint;
%intf2approxs are the approximations for -4G using integers
intf2approxs = Rint./Pint;
%intapproxs are the approximations for G using integers
intapproxs = intf2approxs/(-4);
%error is the difference between intapproxs and G
interror = abs(intapproxs-G);

%We plot intapproxs against G
%Exactly the same as regular approximations (aka without integers)
figure
plot(intapproxs)
hold on
plot(repmat(G,[1,limit]))
title("Type-II integer approximations vs Catalan's constant")
xlabel('n')
hold off

%-----------------------------------------------------------------------
%                   Type-I Hermite-Padé approximations
%-----------------------------------------------------------------------

%Defining the function evaluation f_1(-1)
f1eval =pi/2;

%A will be the vector containing the polynomial evaluations A_n(-1)
A = zeros(1,limit);
%Aint will be the vector containing the polynomial evaluations A_n(-1)
%multiplied by the constant to make them integers
Aint = zeros(1,limit);
%B will be the vector containing the polynomial evaluations B_n(-1)
B = zeros(1,limit);
%Bint will be the vector containing the polynomial evaluations B_n(-1)
%multiplied by the constant to make them integers
Bint = zeros(1,limit);
%C will be the vector containing the polynomial evaluations C_n,n(-1)
C = zeros(1,limit);
%Cint will be the vector containing the polynomial evaluations C_n,n(-1)
%multiplied by the constant to make them integers
Cint = zeros(1,limit);

for n = 1:limit
    noddsl = sym(1:2:2*n-1);
    [A(n),B(n),C(n)] = typeI(n,n);
    Aint(n) = lcm(noddsl)^2*4^(2*n)*A(n);
    Bint(n) = lcm(noddsl)^2*4^(2*n)*B(n);
    Cint(n) = lcm(noddsl)^2*4^(2*n)*C(n);
end

%convergence2 should go to zero; according to our results
convergence2 = abs(A*f1eval - f2eval*B + C);

%We plot convergence2 logarithmically
figure
semilogy(1:limit,convergence2)
title('Logarithmic plot of the Type-I convergence')
xlabel('n')
ylabel('convergence')

%intconvergence2 should NOT go to zero; according to our results
intconvergence2 = abs(Aint*f1eval + f2eval*Bint + Cint);

%We plot intconvergence2 logarithmically
figure
semilogy(1:limit,intconvergence2)
title('Logarithmic plot of the integer Type-I convergence')
xlabel('n')
ylabel('convergence')

%-----------------------------------------------------------------------
%                               Functions
%-----------------------------------------------------------------------

%Calulates polynomial evaluations P_n,m(-1) and R_n,m(-1) for parameters n and m
function [P,R] = typeII(n,m)
    P = 0;
    for k = 0:n+m
        coef = nchoosek(n+m,k)*ratbin(m+k-0.5,m)*ratbin(n+k-0.5,n);
        P = P+coef;
    end
    R = 0;
    for k = 1:n+m
        soem = 0;
        for j = 0:k-1
            term_soem = (-1)^j/(j+0.5)^2;
            soem = soem + term_soem;
        end
        coef = -soem*nchoosek(n+m,k)*ratbin(m+k-0.5,m)*ratbin(n+k-0.5,n);
        R = R+coef;
    end
end

%Calulates polynomial evaluations A_n(-1), B_m(-1) and C_n,m(-1) for parameters n and m
function [A,B,C] = typeI(n,m)
    C = 0;
    ak = zeros(1,n+1);
    bk = zeros(1,m+1);
    for k = 0:m
        a = 0;
        for c = 0:n
            if c ~= k
                coef = nchoosek(n,c)*ratbin(n+0.5+c,n)*nchoosek(m,k)*ratbin(m+n+0.5+k,m)*(-1)^(c+k)*(2*k+1)/(2*(c-k));
                a = a + coef;
            end
        end
        for j = 0:m
            if j ~= k
                coef = nchoosek(n,k)*ratbin(n+0.5+k,n)*nchoosek(m,j)*ratbin(m+n+0.5+j,m)*(-1)^(j+k)*(2*k+1)/(2*(k-j));
                a = a - coef;
            end
        end
        coef = nchoosek(n,k)*ratbin(n+0.5+k,n)*nchoosek(m,k)*ratbin(m+n+0.5+k,m);
        a = a - coef;
        ak(k+1) = a;
        b = nchoosek(n,k)*ratbin(n+0.5+k,n)*nchoosek(m,k)*ratbin(m+n+0.5+k,m)*(2*k+1)/2;
        bk(k+1) = b;
    end
    if m < n
        for k = m+1:n
            a = 0;
            for j = 0:m
                coef = nchoosek(n,k)*ratbin(n+0.5+k,n)*nchoosek(m,j)*ratbin(m+n+0.5+j,m)*(-1)^(k+j)*(2*k+1)/(2*(k-j));
                a = a - coef;
            end
            ak(k+1) = a;
        end
    end
    for k = 0:n
        for l = 0:k-1
            coef = (-1)^l*2*ak(k+1)/(2*k-2*l-1);
            C = C + coef;
        end
    end
    for j = 0:m
        for l = 0:j-1
            coef = (-1)^l*4*bk(j+1)/(2*j-2*l-1)^2;
            C = C + coef; 
        end
    end
    apol = fliplr(ak);
    A = polyval(apol,-1);
    bpol = fliplr(bk);
    B = polyval(bpol,-1);
end

%Calculates nchoosek, but n can be any real number
function bc = ratbin(n,k)
    if n < k
        bc = 0;
    else
        num = 1;
        for s = 0:k-1
            num = num*(n-s);
        end
        bc = num/factorial(k);
    end
end