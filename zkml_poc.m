% PoC - Inferencia verificable con MATLAB
% Álvaro Martínez Alfaro
% Criptografía, Escuela Superior de Ingeniería Informática, UCLM
%
% Nota: Esta PoC implementa la verificabilidad algebraica (R1CS, QAP,
% Schwartz-Zippel) que constituye el núcleo de un SNARK, pero no la
% propiedad de conocimiento cero que lo haría un zk-SNARK. Los PCS, el
% mecanismo de blinding y Fiat-Shamir se desarrollan en la memoria, pero
% su implementación queda fuera del alcance de este desarrollo.

%% 1.1 - Definición de parámetros

p = 65537;
fescala = 16;

%% 1.2 - Definición de la red neuronal

W1 = [0.5, -0.5; 0.375, 0.125];  % capa oculta (2x2)
W2 = [0.5, 0.25];                % capa de salida (1x2)
x  = [0.5; 1.0];                 % vector de entrada

%% 2 - Cuantización de la red neuronal

cuant = @(v) mod(round(v * fescala), p);

x1_c  = cuant(x(1)); x2_c  = cuant(x(2));
fprintf("Entrada cuantizada en F_%d: [%.2f, %.2f] -> [%d, %d].\n", p, x(1), x(2), x1_c, x2_c)

w11_c = cuant(W1(1,1)); w12_q = cuant(W1(1,2));
w21_c = cuant(W1(2,1)); w22_q = cuant(W1(2,2));
fprintf("Capa oculta cuantizada en F_%d: [%.2f, %.2f; %.2f, %.2f] -> [%d, %d; %d, %d].\n", ...
    p, W1(1,1), W1(1,2), W1(2,1), W1(2,2), w11_c, w12_q, w21_c, w22_q)

w31_c = cuant(W2(1)); w32_c = cuant(W2(2));
fprintf("Capa salida cuantizada en F_%d: [%.2f, %.2f] -> [%d, %d].\n", p, W2(1), W2(2), w31_c, w32_c)

%% 3.1 - Inferencia real (punto flotante)

relu = @(v) max(0, v);

s1_float = W1 * x;
fprintf("\ns1_float =\n"); disp(s1_float)

a1_float = relu(s1_float);
fprintf("a1_float =\n"); disp(a1_float)

y_float = W2 * a1_float;
fprintf("y_float = %.10g\n", y_float)

%% 3.2 - Inferencia cuantizada

t1 = mod(w11_c * x1_c, p); % w11 * x1
t2 = mod(w12_q * x2_c, p); % w12 * x2  (negativo en F_p)
t3 = mod(w21_c * x1_c, p); % w21 * x1
t4 = mod(w22_q * x2_c, p); % w22 * x2
fprintf("\nt1=%d  t2=%d  t3=%d  t4=%d\n", t1, t2, t3, t4)

s1_1_q = mod(t1 + t2, p); % pre-activacion neurona 1
s1_2_q = mod(t3 + t4, p); % pre-activacion neurona 2
fprintf("s1_1_q=%d  s1_2_q=%d\n", s1_1_q, s1_2_q)

sgn1 = double(s1_1_q < p/2); % 65473 -> negativo (-64)
sgn2 = double(s1_2_q < p/2); % 80 -> positivo
fprintf("sgn1=%d  sgn2=%d\n", sgn1, sgn2)

a1_1_q = mod(s1_1_q * sgn1, p);
a1_2_q = mod(s1_2_q * sgn2, p);
fprintf("a1_1_q=%d  a1_2_q=%d\n", a1_1_q, a1_2_q)

t5 = mod(w31_c * a1_1_q, p);
t6 = mod(w32_c * a1_2_q, p);
y_q = mod(t5 + t6, p);

y_decuant = y_q / fescala^3;
fprintf("Resultado inferencia real: %g\nResultado inferencia cuantizada: %d (decuantizado: %g)\n", ...
    y_float, y_q, y_decuant)

%% 4 - Construcción del sistema R1CS

% Vector testigo z (20 variables)
z = [1; x1_c; x2_c; y_q; ...
     w11_c; w12_q; w21_c; w22_q; w31_c; w32_c; ...
     t1; t2; t3; t4; ...
     sgn1; a1_1_q; sgn2; a1_2_q; ...
     t5; t6];

fprintf("\nVector testigo: z = ["); fprintf("%d ", z); fprintf("]\n\n")

% Matrices R1CS (m=11 restricciones, n=20 variables)
n = 20; % numero de variables en z
m = 11; % numero de restricciones
A = zeros(m, n); B = zeros(m, n); C = zeros(m, n);

A(1,5)=1;                  B(1,2)=1;    C(1,11)=1; % 1: w11*x1=t1
A(2,6)=1;                  B(2,3)=1;    C(2,12)=1; % 2: w12*x2=t2
A(3,7)=1;                  B(3,2)=1;    C(3,13)=1; % 3: w21*x1=t3
A(4,8)=1;                  B(4,3)=1;    C(4,14)=1; % 4: w22*x2=t4
A(5,11)=1; A(5,12)=1;      B(5,15)=1;   C(5,16)=1; % 5: (t1+t2)*sgn1=a11
A(6,15)=1;                 B(6,15)=1;   C(6,15)=1; % 6: sgn1^2=sgn1
A(7,13)=1; A(7,14)=1;      B(7,17)=1;   C(7,18)=1; % 7: (t3+t4)*sgn2=a12
A(8,17)=1;                 B(8,17)=1;   C(8,17)=1; % 8: sgn2^2=sgn2
A(9,9)=1;                  B(9,16)=1;   C(9,19)=1; % 9: w31*a11=t5
A(10,10)=1;                B(10,18)=1;  C(10,20)=1; % 10: w32*a12=t6
A(11,19)=1; A(11,20)=1;    B(11,1)=1;   C(11,4)=1; % 11: (t5+t6)*1=y_q

% Verificación R1CS
Az = mod(A * z, p);
Bz = mod(B * z, p);
Cz = mod(C * z, p);
izq_r1cs = mod(Az .* Bz, p); % producto de Hadamard
for i = 1:m
    verif = (izq_r1cs(i) == Cz(i));
    marca = 'OK';
    if ~verif, marca = 'ERROR'; end
    fprintf('[%s] Restriccion %2d: %6d x %6d = %6d  (esperado %6d)\n', ...
        marca, i, Az(i), Bz(i), izq_r1cs(i), Cz(i));
end

%% 5 - Transformación a QAP

t_pts = 1:m; % puntos de evaluacion: {1, 2, ..., 11}

A_poli = zeros(1, m);
B_poli = zeros(1, m);
C_poli = zeros(1, m);

for j = 1:n
    if z(j) == 0, continue; end % contribucion nula: z[j]*L = 0

    col_A = A(:, j)';
    col_B = B(:, j)';
    col_C = C(:, j)';

    if any(col_A ~= 0)
        Lj = lagrange_interp_fp(t_pts, col_A, p);
        A_poli = poly_add_fp(A_poli, scalar_mul_poly(z(j), Lj, p), p);
    end
    if any(col_B ~= 0)
        Lj = lagrange_interp_fp(t_pts, col_B, p);
        B_poli = poly_add_fp(B_poli, scalar_mul_poly(z(j), Lj, p), p);
    end
    if any(col_C ~= 0)
        Lj = lagrange_interp_fp(t_pts, col_C, p);
        C_poli = poly_add_fp(C_poli, scalar_mul_poly(z(j), Lj, p), p);
    end
end

% Verificación de la interpolación
fprintf("\n")
for i = 1:m
    Ai = poly_eval_fp(A_poli, i, p);
    Bi = poly_eval_fp(B_poli, i, p);
    Ci = poly_eval_fp(C_poli, i, p);
    verif = (Ai == Az(i)) && (Bi == Bz(i)) && (Ci == Cz(i));
    marca = 'OK';
    if ~verif, marca = 'ERROR'; end
    fprintf('[%s] Punto %2d:  A=%6d (R1CS=%6d)  B=%6d (R1CS=%6d)  C=%6d (R1CS=%6d)\n', ...
        marca, i, Ai, Az(i), Bi, Bz(i), Ci, Cz(i));
end
fprintf('Grado de A(x): %d  B(x): %d  C(x): %d\n', ...
    length(A_poli)-1, length(B_poli)-1, length(C_poli)-1);

% Polinomio de anulación Z(x) = (x-1)(x-2)···(x-11)
Z_poli = [1];
for i = 1:m
    Z_poli = poly_mul_fp(Z_poli, [1, mod(-i, p)], p);
end
fprintf('Grado de Z(x): %d (igual al numero de restricciones)\n\n', length(Z_poli)-1);

% T(x) = A(x)*B(x) - C(x)
AB_poli = poly_mul_fp(A_poli, B_poli, p);       % A(x)*B(x) mod p
neg_C   = scalar_mul_poly(p - 1, C_poli, p);    % -C(x) mod p
T_poli  = poly_add_fp(AB_poli, neg_C, p);        % A(x)*B(x)-C(x) mod p
fprintf('Grado de A(x)*B(x): %d  T(x): %d\n', length(AB_poli)-1, length(T_poli)-1);

% División polinómica: H(x) = T(x) / Z(x)
[H_poli, rem_poli] = poly_div_fp(T_poli, Z_poli, p);
rem_is_zero = isempty(rem_poli) || all(rem_poli == 0);
if rem_is_zero
    fprintf('[OK] T(x) es exactamente divisible por Z(x).\n');
    fprintf('     Grado de H(x) = %d (maximo esperado: %d)\n', ...
        length(H_poli)-1, length(T_poli)-1 - (length(Z_poli)-1));
else
    fprintf('[ERROR] T(x) no es divisible por Z(x). Resto: ');
    fprintf('%d ', rem_poli);
    fprintf('\n');
end

%% 6 - Verificación (Schwartz-Zippel)

r = randi([12, p-1]);
fprintf('\nPunto de evaluacion: r = %d\n', r);

A_r = poly_eval_fp(A_poli, r, p);
B_r = poly_eval_fp(B_poli, r, p);
C_r = poly_eval_fp(C_poli, r, p);
H_r = poly_eval_fp(H_poli, r, p);
Z_r = poly_eval_fp(Z_poli, r, p);

fprintf('A(r) = %d, B(r) = %d, C(r) = %d\n', A_r, B_r, C_r);
fprintf('H(r) = %d, Z(r) = %d\n', H_r, Z_r);

res_izq = mod(mod(A_r * B_r, p) - C_r + p, p);
res_der = mod(H_r * Z_r, p);

fprintf('LHS = A(r)*B(r) - C(r) mod p = %d\n', res_izq);
fprintf('RHS = H(r)*Z(r) mod p = %d\n', res_der);

if res_izq == res_der
    fprintf('[OK] Verificacion exitosa: LHS == RHS\n');
else
    fprintf('[ERROR] Verificacion fallida: LHS != RHS\n');
end

%% Detección de fraude (testigo manipulado)

z_f = z;
z_f(5) = mod(z(5) + 1, p); % modificamos w11 en una unidad

% Recalculamos A(x) con el testigo manipulado
A_poly_f = zeros(1, m);
for j = 1:n
    if z_f(j) == 0, continue; end
    col_A = A(:, j)';
    if any(col_A ~= 0)
        Lj = lagrange_interp_fp(t_pts, col_A, p);
        A_poly_f = poly_add_fp(A_poly_f, scalar_mul_poly(z_f(j), Lj, p), p);
    end
end

A_r_f = poly_eval_fp(A_poly_f, r, p);
lhs_f = mod(mod(A_r_f * B_r, p) - C_r + p, p);
rhs_f = res_der; % la probadora reutiliza H(r)*Z(r) del testigo honesto publicado

if lhs_f ~= rhs_f
    fprintf('\nPrueba falsa DETECTADA en r=%d  (LHS=%d != RHS=%d)\n', r, lhs_f, rhs_f);
else
    fprintf('\n[ADVERTENCIA] Fraude no detectado en r=%d\n', r);
end

% =========================================================================
% FUNCIONES AUXILIARES
% =========================================================================

function inv = modinv_fp(a, p)
% Inverso multiplicativo de a en F_p via teorema de Fermat: a^(p-2) mod p
    inv = powermod(mod(a, p), p-2, p);
end

function r = poly_add_fp(a, b, p)
% Suma de dos polinomios sobre F_p (rellena el mas corto con ceros al frente)
    na = length(a);
    nb = length(b);
    if na >= nb
        r = mod(a + [zeros(1, na - nb), b], p);
    else
        r = mod([zeros(1, nb - na), a] + b, p);
    end
end

function r = poly_mul_fp(a, b, p)
% Producto de polinomios sobre F_p (acumulacion explicita mod p)
    na = length(a); nb = length(b);
    r = zeros(1, na + nb - 1);
    for ii = 1:na
        for jj = 1:nb
            r(ii + jj - 1) = mod(r(ii + jj - 1) + mod(a(ii) * b(jj), p), p);
        end
    end
end

function r = scalar_mul_poly(s, poly, p)
% Multiplica un polinomio por un escalar mod p
    r = mod(mod(s, p) .* poly, p);
end

function [q_out, r_out] = poly_div_fp(num, den, p)
% Division euclidiana de polinomios sobre F_p.
% Implementa el algoritmo de division larga adaptado a aritmetica modular.
% En cada paso, el coeficiente lider del resto se cancela multiplicando el
% divisor por un escalar adecuado y restando.
% En F_p, dividir por el coeficiente lider del divisor equivale a
% multiplicar por su inverso modular.

    % Normalizar: reducir coeficientes mod p y eliminar ceros lideres.
    r_out = mod(num, p);
    while length(r_out) > 1 && r_out(1) == 0
        r_out = r_out(2:end);
    end

    den_c = mod(den, p);
    while length(den_c) > 1 && den_c(1) == 0
        den_c = den_c(2:end);
    end

    % Inverso del coeficiente lider del divisor (se reutiliza en cada paso)
    lc_inv   = modinv_fp(den_c(1), p);
    q_coeffs = zeros(1, max(0, length(r_out) - length(den_c) + 1));
    qi = 0;

    while length(r_out) >= length(den_c)
        % lc = coef_lider(r_out) / coef_lider(den) mod p
        lc = mod(r_out(1) * lc_inv, p);
        qi = qi + 1;
        q_coeffs(qi) = lc;

        % Restar lc * den_c de los primeros coeficientes de r_out
        nd = length(den_c);
        for k = 1:nd
            r_out(k) = mod(r_out(k) - mod(lc * den_c(k), p) + p, p);
        end

        % Eliminar el cero lider resultante de la cancelacion
        while ~isempty(r_out) && r_out(1) == 0
            r_out = r_out(2:end);
        end
    end

    q_out = q_coeffs(1:qi);
    if isempty(q_out), q_out = [0]; end
    % r_out es [] si la division es exacta, o el polinomio resto en otro caso
end

function val = poly_eval_fp(coefs, x, p)
% Evaluacion de polinomio en x sobre F_p usando el metodo de Horner.
% a(1)*x^n + ... + a(n+1) = (...((a(1)*x + a(2))*x + a(3))*x ...) a(n+1))
% En cada paso solo se realiza una multiplicacion y una suma, ambas mod p.
    val = 0;
    x_mod = mod(x, p);
    for k = 1:length(coefs)
        val = mod(val * x_mod + coefs(k), p);
    end
end

function poly = lagrange_interp_fp(t_pts, vals, p)
% Interpolacion de Lagrange sobre F_p.
%
% Dado un conjunto de m puntos distintos t_pts = [t1,...,tm] y sus
% valores vals = [v1,...,vm], construye el unico polinomio de grado
% <= m-1 que pasa por todos ellos.
%
% P(x) = sum_{j=1}^{m} v_j * L_j(x)
% L_j(x) = prod_{k~=j} (x - t_k) / (t_j - t_k)
%
% En F_p, la division por (t_j - t_k) se realiza multiplicando por su
% inverso modular.
    m_loc = length(t_pts);
    poly  = zeros(1, m_loc);

    for j = 1:m_loc
        vj = mod(vals(j), p);
        if vj == 0, continue; end

        % Numerador de L_j(x): prod_{k~=j} (x - t_k)
        num_p = [1];
        for k = 1:m_loc
            if k == j, continue; end
            num_p = poly_mul_fp(num_p, [1, mod(-t_pts(k), p)], p);
        end

        % Denominador: prod_{k~=j} (t_j - t_k) mod p  (escalar)
        denom = 1;
        for k = 1:m_loc
            if k == j, continue; end
            denom = mod(denom * mod(t_pts(j) - t_pts(k), p), p);
        end

        % Escalar el polinomio base por v_j / denom y acumular
        coeff = mod(vj * modinv_fp(denom, p), p);
        poly  = poly_add_fp(poly, scalar_mul_poly(coeff, num_p, p), p);
    end
end
