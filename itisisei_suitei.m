%% 椅子に座った状態で片脚の膝の曲げ伸ばし
clear
close all
clc

%% ==========================
%% パス追加
%% ==========================

addpath('KIT実験')
addpath('富大山内研')

%% ==========================
%% Excel読み込み
%% ==========================

T = readtable('1link_KIT_08.xlsx');

%% 脛部位置
x = T{9:end,26};
y = T{9:end,27};
z = T{9:end,28};

p = [x y z];

%% クォータニオン
xq = T{9:end,22};
yq = T{9:end,23};
zq = T{9:end,24};
wq = T{9:end,25};

q = [wq xq yq zq];

%% 時間
t = T{9:end,2};

%% データ数
N = length(t);

%% ==========================
%% 回転行列へ変換
%% ==========================

R=zeros(3,3,N);

for k=1:N

    quat=quaternion(q(k,:));

    R(:,:,k)=rotmat(quat,'point');

end

%% ==========================
%% 膝中心推定
%% ==========================

pwf0=p(1,:)';
Rwf0=R(:,:,1);

A=[];
B=[];

I=eye(3);

for k=2:N

    pwft=p(k,:)';

    Rwft=R(:,:,k);

    Rtheta=Rwft*Rwf0';

    theta=acos( ...
        max( ...
        min((trace(Rtheta)-1)/2,1), ...
        -1));

    if theta<deg2rad(5)

        continue

    end

    Ak=I-Rtheta;

    Bk=pwft-Rtheta*pwf0;

    A=[A;Ak];

    B=[B;Bk];

end

pwk=lsqminnorm(A,B);

%% ==========================
%% 時間変化する姿勢推定
%% ==========================

Rwk_all=zeros(3,3,N);

for k=2:N

    %% 前フレームとの差分回転

    Rf=R(:,:,k)*R(:,:,k-1)';

    %% 固有値分解

    [V,D]=eig(Rf);

    eigval=diag(D);

    [~,idx]=min(abs(eigval-1));

    nw=real(V(:,idx));

    %% 向き補正

    v=[ ...
        Rf(3,2)-Rf(2,3)
        Rf(1,3)-Rf(3,1)
        Rf(2,1)-Rf(1,2)];

    if dot(nw,v)<0

        nw=-nw;

    end

    ey=nw/norm(nw);

    %% 基準ベクトル

    a=[0;0;1];

    if abs(dot(a,ey))>0.9

        a=[1;0;0];

    end

    ex=a-(a'*ey)*ey;

    ex=ex/norm(ex);

    ez=cross(ex,ey);

    ez=ez/norm(ez);

    Rwk_all(:,:,k)=[ex ey ez];

end

Rwk_all(:,:,1)=Rwk_all(:,:,2);

%% ==========================
%% アニメーション表示
%% ==========================

figure

hold on
grid on
axis equal
view(3)

margin=100;

xlim([min(p(:,1))-margin ...
      max(p(:,1))+margin])

ylim([min(p(:,2))-margin ...
      max(p(:,2))+margin])

zlim([min(p(:,3))-margin ...
      max(p(:,3))+margin])

xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')

title('Time-varying Knee Axis')

%% 軌跡

plot3( ...
    p(:,1),...
    p(:,2),...
    p(:,3),...
    ':k')

%% 膝中心

scatter3( ...
    pwk(1),...
    pwk(2),...
    pwk(3),...
    150,...
    'filled')

L=100;

hx=quiver3(0,0,0,0,0,0,'r','LineWidth',3);
hy=quiver3(0,0,0,0,0,0,'g','LineWidth',3);
hz=quiver3(0,0,0,0,0,0,'b','LineWidth',3);

haxis=plot3(0,0,0,'--g','LineWidth',2);

for k=1:3:N

    Rwk=Rwk_all(:,:,k);

    ex=Rwk(:,1);
    ey=Rwk(:,2);
    ez=Rwk(:,3);

    set(hx,...
        'XData',pwk(1),...
        'YData',pwk(2),...
        'ZData',pwk(3),...
        'UData',L*ex(1),...
        'VData',L*ex(2),...
        'WData',L*ex(3));

    set(hy,...
        'XData',pwk(1),...
        'YData',pwk(2),...
        'ZData',pwk(3),...
        'UData',L*ey(1),...
        'VData',L*ey(2),...
        'WData',L*ey(3));

    set(hz,...
        'XData',pwk(1),...
        'YData',pwk(2),...
        'ZData',pwk(3),...
        'UData',L*ez(1),...
        'VData',L*ez(2),...
        'WData',L*ez(3));

    L2=300;

    set(haxis,...
        'XData',[pwk(1)-L2*ey(1),pwk(1)+L2*ey(1)],...
        'YData',[pwk(2)-L2*ey(2),pwk(2)+L2*ey(2)],...
        'ZData',[pwk(3)-L2*ey(3),pwk(3)+L2*ey(3)])

    drawnow

end