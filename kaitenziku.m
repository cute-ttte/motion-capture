% 椅子に座った状態で片脚の膝の曲げ伸ばし
clear;
close all;
clc;

%% パス追加
addpath('KIT実験');
addpath('富大山内研');

%% Excel読み込み
T = readtable('1link_KIT_08.xlsx');

%% 脛部重心位置
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
    R(:,:,k) = rotmat(quat,'point');
end

%% ==========================
%% 1. 膝関節位置 pwk の推定
%% ==========================
pwf_0 = p(1,:)';
Rwf_0 = R(:,:,1);

A = [];
B = [];
I = eye(3);
theta_array = zeros(N, 1);

for k = 2:N
    pwf_t = p(k,:)';
    Rwf_t = R(:,:,k);

    R_theta = Rwf_t * Rwf_0';
    theta = acos(max(min((trace(R_theta)-1)/2,1),-1));
    theta_array(k) = theta;

    if theta < deg2rad(5)
        continue
    end

    A_k = I - R_theta;
    B_k = pwf_t - R_theta*pwf_0;
    A = [A;A_k];
    B = [B;B_k];
end

pwk = lsqminnorm(A,B);
fprintf('\n推定膝関節位置\n')
fprintf('X = %.3f mm\n',pwk(1))
fprintf('Y = %.3f mm\n',pwk(2))
fprintf('Z = %.3f mm\n',pwk(3))

%% ==========================
%% 2. 膝座標系回転行列 Rwk の推定
%% ==========================
[~, max_idx] = max(theta_array);
R_f = R(:,:,max_idx) * Rwf_0';

[V,D] = eig(R_f);
eigval = diag(D);
[~,idx] = min(abs(eigval-1));
nw = real(V(:,idx));

v_skew = [R_f(3,2)-R_f(2,3); R_f(1,3)-R_f(3,1); R_f(2,1)-R_f(1,2)];
if dot(nw, v_skew) < 0
    nw = -nw;
end

ey = nw/norm(nw);

a = [0;0;1];
if abs(dot(a,ey)) > 0.9
    a = [1;0;0];
end

ex = a - (a'*ey)*ey;
ex = ex/norm(ex);

ez = cross(ex,ey);
ez = ez/norm(ez);

Rwk = [ex ey ez];
disp('推定回転行列 Rwk')
disp(Rwk)

%% ==========================
%% 3. 膝関節角度 θ(t) の算出 (pwk と Rwk の両方を使用)
%% ==========================
theta_signed = zeros(N, 1);
for k = 1:N
    % 1. 世界座標系における現在の脛部マーカー位置
    pwf_t = p(k,:)';
    
    % 2. 膝関節座標系 (ローカル座標系) から見た脛部マーカーの相対位置ベクトルへの変換
    % 剛体変換の逆変換: p_local = R_wk^T * (p_world - p_0)
    p_kf = Rwk' * (pwf_t - pwk);
    
    % 3. Y軸回りの回転とみなして、XZ平面上の座標成分から角度を抽出 (atan2を使用)
    % p_kf(1) がローカルのX座標、p_kf(3) がローカルのZ座標
    x_k = p_kf(1);
    z_k = p_kf(3);
    
    % 角度の算出
    theta_signed(k) = atan2(x_k, z_k);
    
    % 初期姿勢を0度とするためのオフセット調整
    if k == 1
        theta_offset = theta_signed(1);
    end
    theta_signed(k) = theta_signed(k) - theta_offset;
end
% 角度の連続化（±180度の境界をまたぐ場合の不連続性を補正）
theta_signed = unwrap(theta_signed);
%% ==========================
%% 可視化 (3D位置と角度グラフ)
%% ==========================

% --- 3D空間の可視化 ---
figure('Name', '3D Knee Coordinate', 'Position', [100, 100, 600, 500])
hold on; grid on; axis equal; view(3);

plot3(p(:,1), p(:,2), p(:,3), 'LineWidth', 2)
scatter3(pwk(1), pwk(2), pwk(3), 150, 'filled')

L = 100;
quiver3(pwk(1), pwk(2), pwk(3), L*ex(1), L*ex(2), L*ex(3), 'r', 'LineWidth', 2)
quiver3(pwk(1), pwk(2), pwk(3), L*ey(1), L*ey(2), L*ey(3), 'g', 'LineWidth', 2)
quiver3(pwk(1), pwk(2), pwk(3), L*ez(1), L*ez(2), L*ez(3), 'b', 'LineWidth', 2)

xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
legend('Shank trajectory', 'Knee center', 'x-axis', 'y-axis', 'z-axis')
title('Estimated knee coordinate')
rotate3d on

% --- 膝関節角度のグラフ ---
figure('Name', 'Knee Angle', 'Position', [750, 100, 600, 400])
plot(t, rad2deg(theta_signed), 'LineWidth', 2, 'Color', '#D95319')
grid on
xlabel('Time [s]')
ylabel('Flexion Angle [deg]')
title('Knee Joint Angle \theta(t)')
%% ==========================
%% 可視化 (3Dアニメーションと角度グラフのリアルタイム連動)
%% ==========================

% --- グラフ描画の事前準備 ---
% 1. 3D空間アニメーション用のウィンドウ設定
fig3d = figure('Name', '3D Knee Motion Animation', 'Position', [100, 100, 600, 500]);
hold on; grid on; axis equal; view(3);

% 描画範囲（XYZ軸の表示ミニマム・マックス）をデータに合わせて自動調整
margin = 50; % 表示の余白 (mm)
xlim([min([p(:,1); pwk(1)])-margin, max([p(:,1); pwk(1)])+margin]);
ylim([min([p(:,2); pwk(2)])-margin, max([p(:,2); pwk(2)])+margin]);
zlim([min([p(:,3); pwk(3)])-margin, max([p(:,3); pwk(3)])+margin]);

% 【固定オブジェクト】の描画（膝関節の位置・姿勢、および全体の軌跡ガイド線）
h_knee = scatter3(pwk(1), pwk(2), pwk(3), 150, 'filled', 'MarkerFaceColor', 'k'); % 膝中心（黒点）
AxisLength = 100; % 膝座標系の矢印の長さ (mm)
quiver3(pwk(1), pwk(2), pwk(3), AxisLength*ex(1), AxisLength*ex(2), AxisLength*ex(3), 'r', 'LineWidth', 2); % X軸 (赤)
quiver3(pwk(1), pwk(2), pwk(3), AxisLength*ey(1), AxisLength*ey(2), AxisLength*ey(3), 'g', 'LineWidth', 2); % Y軸/回転軸 (緑)
quiver3(pwk(1), pwk(2), pwk(3), AxisLength*ez(1), AxisLength*ez(2), AxisLength*ez(3), 'b', 'LineWidth', 2); % Z軸 (青)
plot3(p(:,1), p(:,2), p(:,3), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1); % 脛部マーカーの全体軌跡（グレー点線）

% 【アニメーション用可変オブジェクト】の初期化 (最初の1フレーム目)
h_shank_point = scatter3(p(1,1), p(1,2), p(1,3), 100, 'filled', 'MarkerFaceColor', '#D95319'); % 現在の脛部マーカー
h_link = plot3([pwk(1), p(1,1)], [pwk(2), p(1,2)], [pwk(3), p(1,3)], 'Color', '#D95319', 'LineWidth', 3); % 膝と脛部を結ぶリンク線
h_time_text = text(pwk(1), pwk(2), max(p(:,3))+30, 'Time: 0.00 s', 'FontSize', 11, 'FontWeight', 'bold'); % 時刻テキスト

xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
legend('Knee Center (pwk)', 'Knee X-axis', 'Knee Y-axis (Rot)', 'Knee Z-axis', 'Shank Trajectory', 'Current Shank Pos', 'Shank Link')
title('3D Knee Joint & Shank Motion')
rotate3d on; % 3Dグラフをマウスで直感的に回転可能にする

% 2. 膝関節角度グラフのウィンドウ設定
fig_ang = figure('Name', 'Knee Angle Plot', 'Position', [750, 100, 600, 400]);
plot(t, rad2deg(theta_signed), 'LineWidth', 2, 'Color', '#0072BD');
grid on; hold on;
xlabel('Time [s]');
ylabel('Flexion Angle [deg]');
title('Knee Joint Angle \theta(t)');
% アニメーションと同期して動く「現在の時刻」を示す赤い縦線
h_xline = xline(t(1), 'r--', 'LineWidth', 1.5);


% --- アニメーションの実行ループ ---
dt = mean(diff(t)); % データの間隔（サンプリング周期）
step = 3;          % 再生スピード・滑らかさの間引き調整（1なら全フレーム描画、数字を大きくすると高速化）

for k = 1:step:N
    % どちらかのウィンドウが閉じられたらループを安全に終了
    if ~ishandle(fig3d) || ~ishandle(fig_ang)
        break;
    end

    % 1. 3D空間プロットの要素を次の時刻のデータに更新
    set(h_shank_point, 'XData', p(k,1), 'YData', p(k,2), 'ZData', p(k,3));
    set(h_link, 'XData', [pwk(1), p(k,1)], 'YData', [pwk(2), p(k,2)], 'ZData', [pwk(3), p(k,3)]);
    set(h_time_text, 'String', sprintf('Time: %.2f s', t(k)));

    % 2. 角度グラフ側の時間縦線を更新
    set(h_xline, 'Value', t(k));

    % 画面の即時書き換えと、実時間に合わせた待機
    drawnow;
    pause(dt * step); 
end