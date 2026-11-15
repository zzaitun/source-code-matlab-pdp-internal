% =========================================================================
% IZRAČUN IN PRIKAZ TABELE GALAT ZA 35 IZBRANIH TOČEK (x, y)
% Model 2D Wave Propagation: ADI vs Implisit Konvensional vs Analitik
% =========================================================================

clear; clc;
fprintf('=========================================================\n');
fprintf(' RINGKASAN GALAT NUMERIK (ADI VS IMPLISIT KONVENSIONAL)\n');
fprintf('=========================================================\n');
% --- 1. Parametri domene in mreže ---
L = 10; M = 10; c = 1.0;
m_mode = 1; n_mode = 1;

Nx = 41; Ny = 41;
dx = L / (Nx - 1); dy = M / (Ny - 1);
x = linspace(0, L, Nx); y = linspace(0, M, Ny);
[X, Y] = meshgrid(x, y);

omega = c * pi * sqrt((m_mode/L)^2 + (n_mode/M)^2);

% --- 2. Časovni parametri ---
dt = 0.08; t_final = 3.0;
Nt = round(t_final / dt);
rx2 = (c * dt / dx)^2; ry2 = (c * dt / dy)^2;

% --- 3. Funkcija analitične rešitve ---
u_exact = @(x_m, y_m, t) sin(m_mode*pi*x_m/L) .* sin(n_mode*pi*y_m/M) .* cos(omega*t);

% Inicializacija matrik
u_imp_nm1 = u_exact(X, Y, 0); u_imp_n = u_exact(X, Y, dt); u_imp_np1 = zeros(Ny, Nx);
u_adi_nm1 = u_exact(X, Y, 0); u_adi_n = u_exact(X, Y, dt); u_adi_np1 = zeros(Ny, Nx);

% --- 4. Časovna iteracija ---
for n = 2:Nt
    % Implisit Konvensional (Proxy)
    for i = 2:Ny-1
        for j = 2:Nx-1
            d2x_n = u_imp_n(i, j+1) - 2*u_imp_n(i, j) + u_imp_n(i, j-1);
            d2y_n = u_imp_n(i+1, j) - 2*u_imp_n(i, j) + u_imp_n(i-1, j);
            u_imp_np1(i, j) = (2*u_imp_n(i, j) - u_imp_nm1(i, j) + 0.5*rx2*d2x_n + 0.5*ry2*d2y_n) / (1 + 0.25*(rx2+ry2));
        end
    end
    u_imp_np1(1,:) = 0; u_imp_np1(end,:) = 0; u_imp_np1(:,1) = 0; u_imp_np1(:,end) = 0;
    u_imp_nm1 = u_imp_n; u_imp_n = u_imp_np1;
    
    % Metode ADI
    u_star = zeros(Ny, Nx);
    A_x = diag((1 + rx2)*ones(Nx-2,1)) + diag(-0.5*rx2*ones(Nx-3,1), 1) + diag(-0.5*rx2*ones(Nx-3,1), -1);
    for i = 2:Ny-1
        d2y_n = u_adi_n(i+1, 2:Nx-1) - 2*u_adi_n(i, 2:Nx-1) + u_adi_n(i-1, 2:Nx-1);
        rhs1 = 2*u_adi_n(i, 2:Nx-1) - u_adi_nm1(i, 2:Nx-1) + ry2*d2y_n;
        u_star(i, 2:Nx-1) = (A_x \ rhs1')';
    end
    
    A_y = diag((1 + ry2)*ones(Ny-2,1)) + diag(-0.5*ry2*ones(Ny-3,1), 1) + diag(-0.5*ry2*ones(Ny-3,1), -1);
    for j = 2:Nx-1
        d2y_nm1 = u_adi_nm1(3:Ny, j) - 2*u_adi_nm1(2:Ny-1, j) + u_adi_nm1(1:Ny-2, j);
        rhs2 = u_star(2:Ny-1, j) - 0.5*ry2*d2y_nm1;
        u_adi_np1(2:Ny-1, j) = A_y \ rhs2;
    end
    u_adi_np1(1,:) = 0; u_adi_np1(end,:) = 0; u_adi_np1(:,1) = 0; u_adi_np1(:,end) = 0;
    u_adi_nm1 = u_adi_n; u_adi_n = u_adi_np1;
end

% --- 5. Izračun napak ---
u_exact_final = u_exact(X, Y, t_final);
err_imp = abs(u_imp_n - u_exact_final);
err_adi = abs(u_adi_n - u_exact_final);

% --- 6. Definiranje 35 izbranih točk (x, y) ---
x_pts = [1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, ...
         2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, ...
         3.75, 3.75, 3.75, ...
         5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, ...
         6.25, 6.25, 6.25, ...
         7.50, 7.50, 7.50, 7.50, 7.50, 7.50, 7.50, ...
         8.75];

y_pts = [1.25, 2.50, 3.75, 5.00, 6.25, 7.50, 8.75, ...
         1.25, 2.50, 3.75, 5.00, 6.25, 7.50, 8.75, ...
         2.50, 5.00, 7.50, ...
         1.25, 2.50, 3.75, 5.00, 6.25, 7.50, 8.75, ...
         2.50, 5.00, 7.50, ...
         1.25, 2.50, 3.75, 5.00, 6.25, 7.50, 8.75, ...
         5.00];

N_points = length(x_pts);
Koordinat_x_y = cell(N_points, 1);
Eksak = zeros(N_points, 1);
Galat_Implisit = zeros(N_points, 1);
Galat_ADI = zeros(N_points, 1);

for k = 1:N_points
    [~, ix] = min(abs(x - x_pts(k)));
    [~, iy] = min(abs(y - y_pts(k)));
    
    Koordinat_x_y{k} = sprintf('(%.2f, %.2f)', x(ix), y(iy));
    Eksak(k) = u_exact_final(iy, ix);
    Galat_Implisit(k) = err_imp(iy, ix);
    Galat_ADI(k) = err_adi(iy, ix);
end

% --- 7. Prikaz tabele ---
Tabel_35_Titik = table(Koordinat_x_y, Eksak, Galat_Implisit, Galat_ADI);
disp(Tabel_35_Titik);