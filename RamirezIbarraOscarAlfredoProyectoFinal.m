%Ramirez Ibarra Oscar Alfredo%
%Proyecto Final

%Filtro analógico de rechazo de banda utilizando 
%una aproximación BButterworth y una estructura Sallen-Key

clear all;
close all;
clc;

%Paso 1. Especificaciones del filtro
Fc1=750;
Fs1=1500;
Fs2=3000;
Fc2=6000;
Ap1=3;
As=60;
Ap2=2.5;

%Paso 2. Normalización a un filtro pasa bajas
Apn=min([Ap1 Ap2]);
Asn=As;

if ((Fs1/Fc1) > (Fc2/Fs2) )
    Fc1_a = Fs1*Fs2/Fc2;
    Fc2_a = Fc2;
else
    Fc1_a = Fc1;
    Fc2_a = Fs1*Fs2/Fc1;
end
Fsn=(Fc2_a - Fc1_a)/(Fs2-Fs1);

%Paso 3. Cálculo del orden del filtro (N)
epsilon=sqrt(10^(Apn/10)-1);
N=log10((10^(Asn/10)-1)/epsilon^2)/(2*log10(Fsn));
N=ceil(N);

%Paso 4. Cálculo de las bicuadráticas
if (mod(N,2))
% N es impar
  Num_N(1,:)=[0 0 epsilon^(-1/N)];
  Den_N(1,:)=[0 1 epsilon^(-1/N)];
  for k=1:(N-1)/2
    Num_N(k+1,:)=[0 0 (epsilon^(-1/N)^2)];
    Den_N(k+1,:)=[1 -2*epsilon^(-1/N)*cos(pi*(N+2*k-1)/(2*N)) (epsilon^(-1/N)^2)];
  end
else
% N es par
  for k=1:N/2
    Num_N(k,:)=[0 0 (epsilon^(-1/N)^2)];
    Den_N(k,:)=[1 -2*epsilon^(-1/N)*cos(pi*(N+2*k-1)/(2*N)) (epsilon^(-1/N)^2)];
  end
end

%Paso 5. Desnormalización del filtro
a = 4*Fc1*Fc2*pi^2;
b = 2*pi*(Fc2 - Fc1);
if (mod(N,2))
    % N es impar
    Num(1,:)=[0 Num_N(1,3)*b 0];
    Den(1,:)=[1 Den_N(1,3)*b a];
    for k=1:(N-1)/2
        p = [1 Num_N(k+1,2)*b (Num_N(k+1,3)*b^2 + 2*a) Num_N(k+1,2)*a*b a^2];
        raices = roots(p);
        Den(2*k,:)=[0 sqrt(Den_N(k+1,3))*b 0];
        Num(2*k,:)=[1 -2*real(raices(1)) abs(raices(1))^2];
        Den(2*k+1,:)=[0 sqrt(Den_N(k+1,3))*b 0];
        Num(2*k+1,:)=[1 -2*real(raices(3)) abs(raices(3))^2];
    end
else
    % N es par
    for k=1:N/2
        p = [1 Num_N(k,2)*b (Num_N(k,3)*b^2 + 2*a) Num_N(k,2)*a*b a^2];
        raices = roots(p);
        Den(2*k-1,:)=[0 sqrt(Den_N(k,3))*b 0];
        Num(2*k-1,:)=[1 -2*real(raices(1)) abs(raices(1))^2];
        Den(2*k,:)=[0 sqrt(Den_N(k,3))*b 0];
        Num(2*k,:)=[1 -2*real(raices(3)) abs(raices(3))^2];
    end
end

% PASO 6. Verificacion de la respuesta en frecuencia
f=logspace(1,5,1000);
s=1j*2*pi*f;

H=1;
  for k=1:N
    H=H.*((Num(k,1)*s.^2 + Num(k,2).*s + Num(k,3))./(Den(k,1)*s.^2 + Den(k,2).*s + Den(k,3)));
  end
figure;
subplot(2,1,1); semilogx(f,20*log10(abs(H))); grid on; hold on;
                %semilogx([10 Fs1 Fs1],[-As -As 0],'r');
                semilogx([Fc1 Fc1 Fc2 Fc2],[-100 -Ap1 -Ap2 -100],'r');
                %semilogx([Fs2 Fs2 10^5],[0 -As -As],'r');
                axis([10 10^5 -100 10]);
subplot(2,1,2); semilogx(f,angle(H)); grid on;

% PASO 7. Determinacion de los componentes de la estructura Tow-Thomas
% Se dijan los valores de C1, C2 y R3
C_ref=10e-9;
R_ref=4700;
b0=Num(k,3);
a1=Den(k,2);
a0=Den(k,3);


for k=1:N
    R(k,3)=R_ref;
    R(k,4)=R_ref;
    C(k,:)=[C_ref C_ref];
    R(k,1)=1/(C(k,1)*a1)
    R(k,2)=1/(a0*R(k,3)*C(k,1)*C(k,2))
    R(k,5)=1/(b0*R(k,3)*C(k,1)*C(k,2))
end







