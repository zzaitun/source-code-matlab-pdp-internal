% =========================================================================
% SIMULASI WAVE PROPAGATION 2D (REFLECTING BOUNDARY)
% Visualisasi 2D Contour dan 3D Surface untuk X-Banner
% =========================================================================

clear; clc; close all;

% Parameter Domain dan Gelombang
Nx = 100; Ny = 100;       % Ukuran grid spasial
c = 1;                    % Kecepatan rambat gelombang
dx = 1; dy = 1;           % Jarak antar grid (spatial step)
dt = 0.5;                 % Langkah waktu (time step)
r = (c * dt / dx)^2;      % Rasio stabilitas

% Inisialisasi Matriks (u_prev: t-1, u: t, u_next: t+1)
u_prev = zeros(Nx, Ny);
u = zeros(Nx, Ny);
u_next = zeros(Nx, Ny);

% Kondisi Awal: Tetesan/Pusat Gelombang di tengah domain (Gaussian Pulse)
[X, Y] = meshgrid(1:Nx, 1:Ny);
x0 = Nx/2; y0 = Ny/2;
u_prev = 0.05 * exp(-((X - x0).^2 + (Y - y0).^2) / 20);
u = u_prev;

% Iterasi Waktu (Menjalankan rambatan gelombang)
% Menggunakan beda hingga standar (proxy visual untuk ADI)
for t = 1:40
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Skema numerik 2D
            u_next(i,j) = 2*u(i,j) - u_prev(i,j) + ...
                          r * (u(i+1,j) - 2*u(i,j) + u(i-1,j) + ...
                               u(i,j+1) - 2*u(i,j) + u(i,j-1));
        end
    end
    % Update waktu
    u_prev = u;
    u = u_next;
end

% ==================== PLOTTING ====================
figure('Position', [100, 100, 1000, 450], 'Color', 'w');

% Grafik 1: 2D Contour (Seperti panel kiri Gambar 2)
subplot(1, 2, 1);
imagesc(u);
colormap(jet);
colorbar;
title('With Reflecting Boundary || 2D Contour', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Distances (x)'); ylabel('Distances (y)');
set(gca, 'YDir', 'normal'); % Membalik sumbu Y agar titik (0,0) di bawah
axis square;

% Grafik 2: 3D Surface (Seperti panel kanan Gambar 2)
subplot(1, 2, 2);
surf(X, Y, u, 'EdgeColor', 'none');
colormap(jet);
title('With Reflecting Boundary || 3D Surface', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('x'); ylabel('y'); zlabel('Amplitudo');
zlim([-0.05 0.05]);
view(-35, 45); % Mengatur sudut pandang kamera 3D
axis square;
grid on;