% 椅子に座った状態で片脚の膝の曲げ伸ばし
clear;
close all;
clc;

%% パス追加
addpath('KIT実験');
addpath('富大山内研');

%% Excel読み込み
T = readtable('1link_KIT_08.xlsx');

%% 脛剛体マーカ重心位置
x = T{9:end,26};
y = T{9:end,27};
z = T{9:end,28};

p = [x y z];

%% クォータニオン
xq = T{9:end,22};
yq = T{9:end,23};
zq = T{9:end,24};
wq = T{9:end,25};

% MATLAB形式 [w x y z]
q = [wq xq yq zq];

%% 時間
t = T{9:end,2};

%% データ数
N = length(t);

%% 回転行列へ変換
R = zeros(3,3,N);

for k = 1:N

    quat = quaternion(q(k,:));

    % 回転行列
    R(:,:,k) = rotmat(quat,'frame');

end

%% 膝関節位置（回転軸）の推定 (最小二乗法)
% 初期状態(t=0)の位置と姿勢
pwf_0 = p(1, :)'; % 3x1ベクトル
Rwf_0 = R(:, :, 1);

% 最小二乗法のための行列AとベクトルBの初期化
% A * pwk = B
A = zeros(3 * N, 3);
B = zeros(3 * N, 1);
I = eye(3);

for k = 1:N
    % 時刻tにおける脛部マーカの位置と姿勢
    pwf_t = p(k, :)';
    Rwf_t = R(:, :, k);
    
    % 相対回転行列 R(θ) の計算 (式 5.3)
    R_theta = Rwf_0' * Rwf_t;
    
    % 方程式の構築: (R(θ) - I) * pwk = R(θ)*pwf(0) - pwf(t)
    A(3*k-2 : 3*k, :) = R_theta - I;
    B(3*k-2 : 3*k, 1) = R_theta * pwf_0 - pwf_t;
end

% 最小二乗法 (バックスラッシュ演算子) で膝関節位置 pwk を一括算出
pwk = A \ B;

%% 結果の計算と表示
disp('=============================================');
fprintf('推定された膝関節(回転軸)の位置 [m]:\n');
fprintf('  x = %.5f\n', pwk(1));
fprintf('  y = %.5f\n', pwk(2));
fprintf('  z = %.5f\n', pwk(3));
disp('=============================================');

%% 座標ごとの時間変化をグラフ表示
figure('Name', 'Position Time Series (Foot vs Knee)', 'Color', 'w');

% グラフのラベルとタイトルの設定
axis_labels = {'X Position [m]', 'Y Position [m]', 'Z Position [m]'};
axis_titles = {'X-axis Time Series', 'Y-axis Time Series', 'Z-axis Time Series'};

for i = 1:3
    subplot(3, 1, i);
    
    % 脛部マーカの各座標の軌跡 (青の実線)
    plot(t, p(:, i), 'b-', 'LineWidth', 1.5); hold on;
    
    % 推定された膝関節位置 (赤の破線、全時刻で一定値)
    plot([t(1), t(end)], [pwk(i), pwk(i)], 'r--', 'LineWidth', 1.5);
    
    % グラフの装飾
    grid on;
    xlabel('Time [s]');
    ylabel(axis_labels{i});
    title(axis_titles{i});
    legend('Foot Marker (p_{wf})', 'Estimated Knee (p_{wk})', 'Location', 'best');
end