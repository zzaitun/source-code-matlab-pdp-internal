% =========================================================================
% SIMULASI WAVE PROPAGATION 2D MENGGUNAKAN METODE ADI
% Kondisi Batas: Absorbing Boundary Condition (Mur's 1st Order ABC)
% Waktu Simulasi: Setelah Gelombang Merambat & Diserap Seluruh Sisi Batas
% =========================================================================

clear; clc; close all;

% --- 1. Parameter Domain dan Fisik ---
L = 10;                    % Panjang domain arah x
M = 10;                    % Panjang domain arah y
Nx = 101;                  % Jumlah titik grid x
Ny = 101;                  % Jumlah titik grid y
dx = L / (Nx - 1);
dy = M / (Ny - 1);

x = linspace(0, L, Nx);
y = linspace(0, M, Ny);
[X, Y] = meshgrid(x, y);

c = 1.0;                   % Kecepatan rambat gelombang
dt = 0.05;                 % Langkah waktu (time step)

% Waktu simulasi disesuaikan (t = 9.53s) seperti pada Gambar 2 proposal
t_final = 25.5;             
Nt = round(t_final / dt);

rx = c * dt / dx;          % Angka Courant arah x
ry = c * dt / dy;          % Angka Courant arah y
rx2 = rx^2;
ry2 = ry^2;

% --- 2. Kondisi Awal (Initial Condition) ---
% Gangguan gelombang Gaussian Pulse dari bagian atas-tengah domain
x0 = 5.0; y0 = 7.0;
u_nm1 = 0.05 * exp(-((X - x0).^2 + (Y - y0).^2) / 0.45);
u_n   = u_nm1;             % Kecepatan awal gelombang = 0
u_np1 = zeros(Ny, Nx);

% --- 3. Penyusunan Matriks Tridiagonal ADI (Interior Grid) ---
% Matriks untuk titik-titik interior (2 sampai N-1)
N_int_x = Nx - 2;
A_x = diag((1 + rx2)*ones(N_int_x, 1)) + ...
      diag(-0.5*rx2*ones(N_int_x - 1, 1), 1) + ...
      diag(-0.5*rx2*ones(N_int_x - 1, 1), -1);

N_int_y = Ny - 2;
A_y = diag((1 + ry2)*ones(N_int_y, 1)) + ...
      diag(-0.5*ry2*ones(N_int_y - 1, 1), 1) + ...
      diag(-0.5*ry2*ones(N_int_y - 1, 1), -1);

% Koefisien Absorbing Boundary Mur Orde-1
coef_x = (rx - 1) / (rx + 1);
coef_y = (ry - 1) / (ry + 1);

% --- 4. Iterasi Waktu Skema ADI ---
for n = 2:Nt
    u_star = zeros(Ny, Nx);
    
    % TAHAP 1: Implisit Arah-x (Untuk titik interior)
    for i = 2:Ny-1
        d2y_n = u_n(i+1, 2:Nx-1) - 2*u_n(i, 2:Nx-1) + u_n(i-1, 2:Nx-1);
        rhs1 = 2*u_n(i, 2:Nx-1) - u_nm1(i, 2:Nx-1) + ry2 * d2y_n;
        u_star(i, 2:Nx-1) = (A_x \ rhs1')';
    end
    
    % TAHAP 2: Implisit Arah-y (Untuk titik interior)
    for j = 2:Nx-1
        d2y_nm1 = u_nm1(3:Ny, j) - 2*u_nm1(2:Ny-1, j) + u_nm1(1:Ny-2, j);
        rhs2 = u_star(2:Ny-1, j) - 0.5 * ry2 * d2y_nm1;
        u_np1(2:Ny-1, j) = A_y \ rhs2;
    end
    
    % --- 5. Penerapan Syarat Batas Serap (Absorbing Boundary Condition) ---
    % Gelombang keluar diserap tanpa dipantulkan kembali ke dalam domain
    % Batas Kiri (x = 0) & Batas Kanan (x = L)
    u_np1(2:Ny-1, 1)  = u_n(2:Ny-1, 2)    + coef_x * (u_np1(2:Ny-1, 2)    - u_n(2:Ny-1, 1));
    u_np1(2:Ny-1, Nx) = u_n(2:Ny-1, Nx-1) + coef_x * (u_np1(2:Ny-1, Nx-1) - u_n(2:Ny-1, Nx));
    
    % Batas Bawah (y = 0) & Batas Atas (y = M)
    u_np1(1, 2:Nx-1)  = u_n(2, 2:Nx-1)    + coef_y * (u_np1(2, 2:Nx-1)    - u_n(1, 2:Nx-1));
    u_np1(Ny, 2:Nx-1) = u_n(Ny-1, 2:Nx-1) + coef_y * (u_np1(Ny-1, 2:Nx-1) - u_n(Ny, 2:Nx-1));
    
    % Penanganan 4 Titik Sudut Domain
    u_np1(1, 1)   = 0.5 * (u_np1(1, 2) + u_np1(2, 1));
    u_np1(1, Nx)  = 0.5 * (u_np1(1, Nx-1) + u_np1(2, Nx));
    u_np1(Ny, 1)  = 0.5 * (u_np1(Ny, 2) + u_np1(Ny-1, 1));
    u_np1(Ny, Nx) = 0.5 * (u_np1(Ny, Nx-1) + u_np1(Ny-1, Nx));
    
    % Update tingkat waktu
    u_nm1 = u_n;
    u_n   = u_np1;
end

% --- 6. Visualisasi Grafik 2D Contour dan 3D Surface ---
figure('Position', [100, 100, 1100, 460], 'Color', 'w');

% Panel Kiri: 2D Contour Plot
subplot(1, 2, 1);
imagesc(x, y, u_n);
colormap(jet);
colorbar;
caxis([-0.01 0.01]); % Rentang skala warna sesuai Gambar 2 proposal
title(sprintf('With Absorbing Boundary || Time = %.2f', t_final), ...
      'FontWeight', 'bold', 'FontSize', 11);
xlabel('x'); ylabel('y');
set(gca, 'YDir', 'normal');
axis square;

% Panel Kanan: 3D Surface Plot
subplot(1, 2, 2);
surf(X, Y, u_n, 'EdgeColor', 'k', 'EdgeAlpha', 0.15);
colormap(jet);
title(sprintf('With Absorbing Boundary || Time = %.2f', t_final), ...
      'FontWeight', 'bold', 'FontSize', 11);
xlabel('x'); ylabel('y'); zlabel('Amplitudo');
zlim([-0.01 0.04]);  % Rentang sumbu Z sesuai Gambar 2 proposal
caxis([-0.01 0.01]);
view(-38, 42);       % Sudut pandang kamera 3D
axis square;
grid on;