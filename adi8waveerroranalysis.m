% =========================================================================
% PERBANDINGAN RMSE & EFISIENSI KOMPUTASI: IMPLISIT KONVENSIONAL VS ADI
% Variasi Pembagian Grid: N = 10, 20, 40, 50, 80, 100
% Fungsi Analitik: f(x,y,t) = A * sin(kx*x + ky*y - omega*t)
% =========================================================================

clear; clc; close all;

% --- 1. Parameter Fisik & Gelombang ---
L = 10;                % Panjang domain x
M = 10;                % Panjang domain y
A = 1.0;               % Amplitudo gelombang
c = 1.0;               % Kecepatan gelombang
kx = 2 * pi / L;       % Bilangan gelombang arah x
ky = 2 * pi / M;       % Bilangan gelombang arah y
omega = c * sqrt(kx^2 + ky^2); % Frekuensi sudut

% Fungsi Analitik Eksak
f_exact = @(x, y, t) A * sin(kx * x + ky * y - omega * t);

% --- 2. Parameter Simulasi ---
grids = [10, 20, 40, 50, 80, 100]; % Variasi pembagian grid
num_grids = length(grids);
dt = 0.02;                         % Langkah waktu (time step)
t_final = 2.0;                     % Waktu akhir simulasi
Nt = round(t_final / dt);          % Total iterasi waktu

% Variabel penyimpan hasil
RMSE_Konvensional = zeros(num_grids, 1);
RMSE_ADI          = zeros(num_grids, 1);
Waktu_Konvensional= zeros(num_grids, 1);
Waktu_ADI         = zeros(num_grids, 1);

% --- 3. Loop Pengujian untuk Setiap Ukuran Grid ---
for g = 1:num_grids
    N = grids(g);
    Nx = N + 1;
    Ny = N + 1;
    dx = L / N;
    dy = M / N;
    
    x = linspace(0, L, Nx);
    y = linspace(0, M, Ny);
    [X, Y] = meshgrid(x, y);
    
    rx2 = (c * dt / dx)^2;
    ry2 = (c * dt / dy)^2;
    
    % =====================================================================
    % A. METODE IMPLISIT KONVENSIONAL (FULL IMPLICIT 2D MATRIX)
    % Skema: (I - rx2*Dxx - ry2*Dyy) u^{n+1} = 2*u^n - u^{n-1}
    % =====================================================================
    tic;
    nx_i = Nx - 2;
    ny_i = Ny - 2;
    
    % Konstruksi Matriks Sistem Besar 2D menggunakan Produk Kronecker
    Ix = speye(nx_i);
    Iy = speye(ny_i);
    ex = ones(nx_i, 1);
    ey = ones(ny_i, 1);
    Dxx = spdiags([ex, -2*ex, ex], [-1, 0, 1], nx_i, nx_i);
    Dyy = spdiags([ey, -2*ey, ey], [-1, 0, 1], ny_i, ny_i);
    
    M_imp = speye(nx_i * ny_i) - rx2 * kron(Ix, Dyy) - ry2 * kron(Dxx, Iy);
    [L_lu, U_lu, P_lu, Q_lu] = lu(M_imp); % Dekomposisi LU Sparse
    
    u_imp_nm1 = f_exact(X, Y, 0);
    u_imp_n   = f_exact(X, Y, dt);
    
    for n = 2:Nt
        t_next = n * dt;
        u_ex_next = f_exact(X, Y, t_next);
        
        rhs = 2 * u_imp_n(2:Ny-1, 2:Nx-1) - u_imp_nm1(2:Ny-1, 2:Nx-1);
        
        % Kontribusi syarat batas Dirichlet bergantung waktu
        rhs(1, :)   = rhs(1, :)   + ry2 * u_ex_next(1, 2:Nx-1);
        rhs(end, :) = rhs(end, :) + ry2 * u_ex_next(end, 2:Nx-1);
        rhs(:, 1)   = rhs(:, 1)   + rx2 * u_ex_next(2:Ny-1, 1);
        rhs(:, end) = rhs(:, end) + rx2 * u_ex_next(2:Ny-1, end);
        
        % Selesaikan sistem linear besar
        sol_vec = Q_lu * (U_lu \ (L_lu \ (P_lu * rhs(:))));
        
        u_imp_np1 = u_ex_next;
        u_imp_np1(2:Ny-1, 2:Nx-1) = reshape(sol_vec, [ny_i, nx_i]);
        
        u_imp_nm1 = u_imp_n;
        u_imp_n   = u_imp_np1;
    end
    Waktu_Konvensional(g) = toc;
    
    % =====================================================================
    % B. METODE ADI (ALTERNATING DIRECTION IMPLICIT ORDE-2)
    % Skema Faktorisasi Konsisten:
    % (1 - rx2/2 d2x)(1 - ry2/2 d2y) u^{n+1} = 2u^n - (1 - rx2/2 d2x)(1 - ry2/2 d2y) u^{n-1}
    % =====================================================================
    tic;
    % Matriks Tridiagonal 1D (Sangat Ringan)
    Ax = diag((1 + rx2)*ones(nx_i, 1)) + diag(-0.5*rx2*ones(nx_i-1, 1), 1) + diag(-0.5*rx2*ones(nx_i-1, 1), -1);
    Ay = diag((1 + ry2)*ones(ny_i, 1)) + diag(-0.5*ry2*ones(ny_i-1, 1), 1) + diag(-0.5*ry2*ones(ny_i-1, 1), -1);
    
    u_adi_nm1 = f_exact(X, Y, 0);
    u_adi_n   = f_exact(X, Y, dt);
    
    for n = 2:Nt
        t_next = n * dt;
        u_ex_next = f_exact(X, Y, t_next);
        
        % Hitung operator kanan pada tingkat waktu n-1
        w = u_adi_nm1;
        w(2:Ny-1, :) = u_adi_nm1(2:Ny-1, :) - 0.5*ry2 * (u_adi_nm1(3:Ny, :) - 2*u_adi_nm1(2:Ny-1, :) + u_adi_nm1(1:Ny-2, :));
        
        z = w;
        z(:, 2:Nx-1) = w(:, 2:Nx-1) - 0.5*rx2 * (w(:, 3:Nx) - 2*w(:, 2:Nx-1) + w(:, 1:Nx-2));
        
        % Syarat batas variabel perantara u_star = (1 - ry2/2 d2y) u^{n+1}
        d2y_ex = f_exact(X, Y+dy, t_next) - 2*f_exact(X, Y, t_next) + f_exact(X, Y-dy, t_next);
        u_star_bc = u_ex_next - 0.5 * ry2 * d2y_ex;
        u_star = u_star_bc;
        
        % Tahap 1: Sapuan Horizontal (Arah-x)
        for i = 2:Ny-1
            rhs1 = 2 * u_adi_n(i, 2:Nx-1) - z(i, 2:Nx-1);
            rhs1(1)   = rhs1(1)   + 0.5 * rx2 * u_star_bc(i, 1);
            rhs1(end) = rhs1(end) + 0.5 * rx2 * u_star_bc(i, end);
            u_star(i, 2:Nx-1) = (Ax \ rhs1')';
        end
        
        % Tahap 2: Sapuan Vertikal (Arah-y)
        u_adi_np1 = u_ex_next;
        for j = 2:Nx-1
            rhs2 = u_star(2:Ny-1, j);
            rhs2(1)   = rhs2(1)   + 0.5 * ry2 * u_ex_next(1, j);
            rhs2(end) = rhs2(end) + 0.5 * ry2 * u_ex_next(end, j);
            u_adi_np1(2:Ny-1, j) = Ay \ rhs2;
        end
        
        u_adi_nm1 = u_adi_n;
        u_adi_n   = u_adi_np1;
    end
    Waktu_ADI(g) = toc;
    
    % =====================================================================
    % C. PERHITUNGAN RMSE TERHADAP FUNGSI ANALITIK
    % =====================================================================
    u_analitik = f_exact(X, Y, Nt * dt);
    
    RMSE_Konvensional(g) = sqrt(mean((u_imp_n(:) - u_analitik(:)).^2));
    RMSE_ADI(g)          = sqrt(mean((u_adi_n(:) - u_analitik(:)).^2));
end

% --- 4. Tampilkan Tabel Hasil di Command Window ---
Grid = grids';
Ratio_of_acuration = RMSE_Konvensional ./ RMSE_ADI;

Tabel_RMSE = table(Grid, RMSE_Konvensional, RMSE_ADI, Ratio_of_acuration, Waktu_Konvensional, Waktu_ADI);
disp('=========================================================================================');
disp('   TABEL PERBANDINGAN RMSE & WAKTU KOMPUTASI (IMPLISIT KONVENSIONAL VS METODE ADI)');
disp('=========================================================================================');
disp(Tabel_RMSE);

% --- 5. Visualisasi Grafik Perbandingan RMSE & Waktu Komputasi ---
figure('Position', [100, 100, 1100, 460], 'Color', 'w');

% Panel Kiri: Grafik Konvergensi RMSE (Skala Logaritmik)
subplot(1, 2, 1);
semilogy(grids, RMSE_Konvensional, 'r--o', 'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'DisplayName', 'Implisit Konvensional'); hold on;
semilogy(grids, RMSE_ADI, 'b-s', 'LineWidth', 2.0, 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'DisplayName', 'Metode ADI (Orde-2)');
title(' RMSE vs Grid Dimension', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Pembagian Grid (N \times N)');
ylabel('RMSE (Skala Logaritmik)');
xticks(grids);
grid on;
legend('Location', 'southwest');

% Panel Kanan: Grafik Perbandingan Waktu Komputasi (CPU Time)
subplot(1, 2, 2);
plot(grids, Waktu_Konvensional, 'r--o', 'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'DisplayName', 'Implisit Konvensional (Matriks 2D)'); hold on;
plot(grids, Waktu_ADI, 'b-s', 'LineWidth', 2.0, 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'DisplayName', 'Metode ADI (Tridiagonal 1D)');
title('Time Computation vs grid dimension', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Grid (N x N)');
ylabel('Computation (Second)');
xticks(grids);
grid on;
legend('Location', 'northwest');