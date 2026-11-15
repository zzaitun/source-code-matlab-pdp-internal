% =========================================================================
% ANALISIS GALAT MODUL 2D WAVE PROPAGATION
% Perbandingan: Beda Hingga Konvensional (Eksplisit) vs ADI vs Solusi Analitik
% Irisan: x = L/2 (Vertikal) dan y = M/2 (Horizontal)
% =========================================================================

clear; clc; close all;

% --- 1. Parameter Fisik dan Grid Domain ---
L = 10;            % Panjang domain x
M = 10;            % Panjang domain y
c = 1.0;           % Kecepatan gelombang
m_mode = 1;        % Moda gelombang x
n_mode = 1;        % Moda gelombang y

Nx = 51;           % Jumlah titik grid arah x (Ganjil agar ada titik tengah tepat)
Ny = 51;           % Jumlah titik grid arah y
dx = L / (Nx - 1);
dy = M / (Ny - 1);

x = linspace(0, L, Nx);
y = linspace(0, M, Ny);
[X, Y] = meshgrid(x, y);

% Frekuensi alami (omega) solusi analitik
omega = c * pi * sqrt((m_mode/L)^2 + (n_mode/M)^2);

% --- 2. Parameter Waktu ---
dt = 0.04;                  % Langkah waktu
t_final = 3.0;              % Waktu akhir simulasi
Nt = round(t_final / dt);   % Jumlah iterasi waktu
rx2 = (c * dt / dx)^2;       % Rasio Courant kuadrat (x)
ry2 = (c * dt / dy)^2;       % Rasio Courant kuadrat (y)

% --- 3. Inisialisasi Matriks Solusi ---
% Solusi Eksak
u_exact = @(x_m, y_m, t) sin(m_mode*pi*x_m/L) .* sin(n_mode*pi*y_m/M) .* cos(omega*t);

% Beda Hingga Konvensional (Eksplisit)
u_exp_nm1 = u_exact(X, Y, 0);
u_exp_n   = u_exact(X, Y, dt);
u_exp_np1 = zeros(Ny, Nx);

% Metode ADI (3-Tingkat Waktu Fairweather-Mitchell)
u_adi_nm1 = u_exact(X, Y, 0);
u_adi_n   = u_exact(X, Y, dt);
u_adi_np1 = zeros(Ny, Nx);

% --- 4. Iterasi Waktu (Time Stepping) ---
for n = 2:Nt
    t_curr = n * dt;
    
    % A. BEDA HINGGA KONVENSIONAL (EKSPLISIT 2D)
    for i = 2:Ny-1
        for j = 2:Nx-1
            d2x = u_exp_n(i, j+1) - 2*u_exp_n(i, j) + u_exp_n(i, j-1);
            d2y = u_exp_n(i+1, j) - 2*u_exp_n(i, j) + u_exp_n(i-1, j);
            u_exp_np1(i, j) = 2*u_exp_n(i, j) - u_exp_nm1(i, j) + rx2*d2x + ry2*d2y;
        end
    end
    % Syarat batas Dirichlet
    u_exp_np1(1,:) = 0; u_exp_np1(end,:) = 0;
    u_exp_np1(:,1) = 0; u_exp_np1(:,end) = 0;
    
    % Update Konvensional
    u_exp_nm1 = u_exp_n;
    u_exp_n   = u_exp_np1;
    
    % B. METODE ADI (ALTERNATING DIRECTION IMPLICIT)
    % Tahap 1: Implisit Arah-x (Mencari u_star)
    u_star = zeros(Ny, Nx);
    A_x = diag((1 + rx2)*ones(Nx-2,1)) + diag(-0.5*rx2*ones(Nx-3,1), 1) + diag(-0.5*rx2*ones(Nx-3,1), -1);
    
    for i = 2:Ny-1
        d2y_n = u_adi_n(i+1, 2:Nx-1) - 2*u_adi_n(i, 2:Nx-1) + u_adi_n(i-1, 2:Nx-1);
        rhs1 = 2*u_adi_n(i, 2:Nx-1) - u_adi_nm1(i, 2:Nx-1) + ry2*d2y_n;
        u_star(i, 2:Nx-1) = (A_x \ rhs1')';
    end
    
    % Tahap 2: Implisit Arah-y (Mencari u^{n+1})
    A_y = diag((1 + ry2)*ones(Ny-2,1)) + diag(-0.5*ry2*ones(Ny-3,1), 1) + diag(-0.5*ry2*ones(Ny-3,1), -1);
    
    for j = 2:Nx-1
        d2y_nm1 = u_adi_nm1(3:Ny, j) - 2*u_adi_nm1(2:Ny-1, j) + u_adi_nm1(1:Ny-2, j);
        rhs2 = u_star(2:Ny-1, j) - 0.5*ry2*d2y_nm1;
        u_adi_np1(2:Ny-1, j) = A_y \ rhs2;
    end
    
    % Syarat batas Dirichlet
    u_adi_np1(1,:) = 0; u_adi_np1(end,:) = 0;
    u_adi_np1(:,1) = 0; u_adi_np1(:,end) = 0;
    
    % Update ADI
    u_adi_nm1 = u_adi_n;
    u_adi_n   = u_adi_np1;
end

% --- 5. Evaluasi Solusi Analitik & Galat Mutlak ---
u_exact_final = u_exact(X, Y, t_final);

err_exp = abs(u_exp_n - u_exact_final);
err_adi = abs(u_adi_n - u_exact_final);

% Indeks Titik Tengah
mid_x = (Nx + 1) / 2; % Titik x = L/2
mid_y = (Ny + 1) / 2; % Titik y = M/2

% --- 6. Plotting Irisan dan Galat ---
figure('Position', [100, 80, 1100, 750], 'Color', 'w');

% Plot 1: Irisan Vertikal x = L/2 (Profil Amplitudo)
subplot(2, 2, 1);
plot(y, u_exact_final(:, mid_x), 'k-', 'LineWidth', 2, 'DisplayName', 'Analitik'); hold on;
plot(y, u_exp_n(:, mid_x), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Konvensional (Eksplisit)');
plot(y, u_adi_n(:, mid_x), 'b-.', 'LineWidth', 1.5, 'DisplayName', 'Metode ADI');
title(['Profil Amplitudo Gelombang pada x = L/2 (t = ', num2str(t_final), 's)'], 'FontWeight', 'bold');
xlabel('Posisi y'); ylabel('Amplitudo u');
grid on; legend('Location', 'best');

% Plot 2: Galat Mutlak Irisan Vertikal x = L/2
subplot(2, 2, 2);
plot(y, err_exp(:, mid_x), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Galat Konvensional'); hold on;
plot(y, err_adi(:, mid_x), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Galat ADI');
title('Galat Mutlak |u_{num} - u_{eksak}| pada x = L/2', 'FontWeight', 'bold');
xlabel('Posisi y'); ylabel('Galat Mutlak');
grid on; legend('Location', 'best');

% Plot 3: Irisan Horizontal y = M/2 (Profil Amplitudo)
subplot(2, 2, 3);
plot(x, u_exact_final(mid_y, :), 'k-', 'LineWidth', 2, 'DisplayName', 'Analitik'); hold on;
plot(x, u_exp_n(mid_y, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Konvensional (Eksplisit)');
plot(x, u_adi_n(mid_y, :), 'b-.', 'LineWidth', 1.5, 'DisplayName', 'Metode ADI');
title(['Profil Amplitudo Gelombang pada y = M/2 (t = ', num2str(t_final), 's)'], 'FontWeight', 'bold');
xlabel('Posisi x'); ylabel('Amplitudo u');
grid on; legend('Location', 'best');

% Plot 4: Galat Mutlak Irisan Horizontal y = M/2
subplot(2, 2, 4);
plot(x, err_exp(mid_y, :), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Galat Konvensional'); hold on;
plot(x, err_adi(mid_y, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Galat ADI');
title('Galat Mutlak |u_{num} - u_{eksak}| pada y = M/2', 'FontWeight', 'bold');
xlabel('Posisi x'); ylabel('Galat Mutlak');
grid on; legend('Location', 'best');

% --- 7. Output Ringkasan Statistik Galat di Command Window ---
fprintf('=========================================================\n');
fprintf(' RINGKASAN GALAT NUMERIK PADA IRISAN GARIS TENGAH\n');
fprintf('=========================================================\n');
fprintf('Irisan x = L/2 (Garis Tengah Vertikal):\n');
fprintf('  - Max Galat Konvensional : %.6e\n', max(err_exp(:, mid_x)));
fprintf('  - Max Galat ADI          : %.6e\n', max(err_adi(:, mid_x)));
fprintf('---------------------------------------------------------\n');
fprintf('Irisan y = M/2 (Garis Tengah Horizontal):\n');
fprintf('  - Max Galat Konvensional : %.6e\n', max(err_exp(mid_y, :)));
fprintf('  - Max Galat ADI          : %.6e\n', max(err_adi(mid_y, :)));
fprintf('=========================================================\n');