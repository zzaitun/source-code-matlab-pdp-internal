% =========================================================================
% SIMULASI WAVE PROPAGATION 2D MENGGUNAKAN METODE ADI
% Kondisi Batas: Reflecting Boundary (Neumann Boundary Condition: du/dn = 0)
% Waktu Simulasi: Setelah Gelombang Terefleksi dari Seluruh Sisi Batas
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

% Waktu simulasi diatur cukup lama (t = 8.5s) agar gelombang 
% sempat menabrak keempat dinding batas dan terefleksi penuh ke dalam
t_final = 25.5;             
Nt = round(t_final / dt);

rx2 = (c * dt / dx)^2;     % Rasio Courant kuadrat arah x
ry2 = (c * dt / dy)^2;     % Rasio Courant kuadrat arah y

% --- 2. Kondisi Awal (Initial Condition) ---
% Gangguan gelombang berbentuk Gaussian Pulse sedikit digeser dari pusat 
% agar refleksi dari setiap dinding terlihat dinamis
x0 = 4.5; y0 = 6.5;
u_nm1 = 0.05 * exp(-((X - x0).^2 + (Y - y0).^2) / 0.4);
u_n   = u_nm1;             % Kecepatan awal gelombang = 0 (du/dt = 0)
u_np1 = zeros(Ny, Nx);

% --- 3. Penyusunan Matriks Tridiagonal ADI (Dengan Reflecting Boundary) ---
% Untuk syarat batas reflektif (Neumann: du/dn = 0), baris pertama dan 
% terakhir matriks diubah menggunakan pendekatan ghost-point orde-2.

% Matriks A_x (Implisit Arah-x)
main_x = (1 + rx2) * ones(Nx, 1);
sub_x  = -0.5 * rx2 * ones(Nx - 1, 1);
sup_x  = -0.5 * rx2 * ones(Nx - 1, 1);
sup_x(1)   = -rx2;         % Refleksi batas kiri (x = 0)
sub_x(end) = -rx2;         % Refleksi batas kanan (x = L)
A_x = diag(main_x) + diag(sup_x, 1) + diag(sub_x, -1);

% Matriks A_y (Implisit Arah-y)
main_y = (1 + ry2) * ones(Ny, 1);
sub_y  = -0.5 * ry2 * ones(Ny - 1, 1);
sup_y  = -0.5 * ry2 * ones(Ny - 1, 1);
sup_y(1)   = -ry2;         % Refleksi batas bawah (y = 0)
sub_y(end) = -ry2;         % Refleksi batas atas (y = M)
A_y = diag(main_y) + diag(sup_y, 1) + diag(sub_y, -1);

% --- 4. Iterasi Waktu Skema ADI ---
for n = 2:Nt
    u_star = zeros(Ny, Nx);
    
    % TAHAP 1: Implisit Arah-x (Mencari variabel perantara u_star)
    for i = 1:Ny
        % Hitung turunan kedua arah-y dengan syarat batas reflektif
        if i == 1
            d2y_n = 2 * (u_n(2, :) - u_n(1, :));
        elseif i == Ny
            d2y_n = 2 * (u_n(Ny-1, :) - u_n(Ny, :));
        else
            d2y_n = u_n(i+1, :) - 2*u_n(i, :) + u_n(i-1, :);
        end
        
        rhs1 = 2*u_n(i, :) - u_nm1(i, :) + ry2 * d2y_n;
        u_star(i, :) = (A_x \ rhs1')';
    end
    
    % TAHAP 2: Implisit Arah-y (Mencari solusi u^{n+1})
    for j = 1:Nx
        % Hitung turunan kedua arah-y pada tingkat waktu n-1 (reflektif)
        d2y_nm1 = zeros(Ny, 1);
        d2y_nm1(1)      = 2 * (u_nm1(2, j) - u_nm1(1, j));
        d2y_nm1(2:Ny-1) = u_nm1(3:Ny, j) - 2*u_nm1(2:Ny-1, j) + u_nm1(1:Ny-2, j);
        d2y_nm1(Ny)     = 2 * (u_nm1(Ny-1, j) - u_nm1(Ny, j));
        
        rhs2 = u_star(:, j) - 0.5 * ry2 * d2y_nm1;
        u_np1(:, j) = A_y \ rhs2;
    end
    
    % Update tingkat waktu
    u_nm1 = u_n;
    u_n   = u_np1;
end

% --- 5. Visualisasi Grafik 2D Contour dan 3D Surface ---
figure('Position', [100, 100, 1100, 460], 'Color', 'w');

% Panel Kiri: 2D Contour Plot
subplot(1, 2, 1);
imagesc(x, y, u_n);
colormap(jet);
colorbar;
caxis([-0.015 0.015]); % Mengatur rentang skala warna agar pola interferensi jelas
title(sprintf('With Reflecting Boundary || Time = %.2f', t_final), ...
      'FontWeight', 'bold', 'FontSize', 11);
xlabel('x'); ylabel('y');
set(gca, 'YDir', 'normal');
axis square;

% Panel Kanan: 3D Surface Plot
subplot(1, 2, 2);
surf(X, Y, u_n, 'EdgeColor', 'k', 'EdgeAlpha', 0.15);
colormap(jet);
title(sprintf('With Reflecting Boundary || Time = %.2f', t_final), ...
      'FontWeight', 'bold', 'FontSize', 11);
xlabel('x'); ylabel('y'); zlabel('Amplitudo');
zlim([-0.03 0.03]);
caxis([-0.015 0.015]);
view(-38, 42); % Sudut pandang kamera 3D
axis square;
grid on;