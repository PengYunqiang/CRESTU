# CRESTU 理论与技术参考手册

**适用源码基线：** `Source Code/` 下 55 个 MATLAB 文件（2026-08-25 规范化版本）  
**相位约定：** $\exp(\mathrm{i}\omega t)$  
**单位制：** SI  
**主要作者：** Yunqiang Peng，Zhentao Jiang（SJTU）

---

## 0. 手册范围、核查方法与实现边界

本手册以源码为唯一实现判据。每一条标为“源码实现式”的关键公式之后，均紧跟实际 MATLAB 核心切片；代码块上方给出文件与冻结行号。理论上应成立、但当前目录没有对应求解路径的关系，明确标为“理论完备式/当前未闭合”，不得据此误判为已实现功能。

源码审计得到四项重要结论：

1. 主求解器采用以边界势 $\Phi$ 为未知量的**直接 Rankine 边界积分方程**；经典 Hess–Smith 间接源法的未知量是源强 $\sigma$。两者共享 Rankine 核与面元解析积分，但未知量、矩阵物理意义不同。
2. `hess_smith_panel_velocity.m` 提供单位常源面元诱导速度解析核，但当前 `assemble_rankine_matrix.m` 实际调用的是 `rankine_panel_integrals.m`。
3. `compute_wave_excitation.m` 能由给定的入射势和绕射势积分一阶载荷；本目录未提供把 $-\partial_n\Phi_0$ 组装成绕射方程右端并求解 $\Phi_7$ 的独立入口。
4. `compute_drift_nearfield.m` 将平均漂移力拆成四个可审计数组，但其“平移—梯度项”和“转动项”是裁剪实现；完整 Pinkster 压力梯度因子及惯性力矩耦合未在该文件中出现。

---

## 1. 坐标系、自由度、常数与符号

### 1.1 全局与局部坐标系

全局大地坐标系记为 $O-xyz$：$x$、$y$ 位于静水面内，$z$ 轴竖直向上，静水面为 $z=0$；有限水深海底为 $z=-h$。第 $k$ 个浮体的局部坐标系记为 $O_k-x_ky_kz_k$，其参考位置（通常为重心）是 $\boldsymbol{x}_{G,k}$。若局部网格绕 $z$ 轴转角为 $\psi_k$，则

\[
\boldsymbol{x}=\boldsymbol{x}_{0,k}+\mathbf R_z(\psi_k)\boldsymbol{x}_k,
\qquad
\mathbf R_z(\psi_k)=
\begin{bmatrix}
\cos\psi_k&-\sin\psi_k&0\\
\sin\psi_k& \cos\psi_k&0\\
0&0&1
\end{bmatrix}. \tag{1-1}
\]

【源码对照：`Source Code/2.Mesh/transform_body_mesh.m`，第 29–36 行】

```matlab
    angle = body_cfg.heading_deg * pi / 180;
    Rz = [cos(angle), -sin(angle), 0; sin(angle), cos(angle), 0; 0, 0, 1];
    x0 = reshape(body_cfg.position, 1, 3);
    mesh_global.vertices = mesh_local.vertices;
    for p = 1:mesh_local.n_panels
        mesh_global.vertices(p, :, :) = squeeze(mesh_local.vertices(p, :, :)) * Rz' + x0;
    end
    mesh_global.centers = mesh_local.centers * Rz' + x0;
```

### 1.2 六自由度与广义法向量

第 $k$ 个浮体的复运动幅值按

\[
\boldsymbol{\xi}_k=
[\xi_1,\xi_2,\xi_3,\xi_4,\xi_5,\xi_6]^{\mathsf T}
=[\text{surge},\text{sway},\text{heave},\text{roll},\text{pitch},\text{yaw}]^{\mathsf T}. \tag{1-2}
\]

面元中心相对重心的位置是 $\boldsymbol r=\boldsymbol x-\boldsymbol x_G$。微小刚体位移为

\[
\delta\boldsymbol x
=\boldsymbol\xi_T+\boldsymbol\xi_R\times\boldsymbol r,
\qquad
\delta\boldsymbol x\cdot\boldsymbol n
=\boldsymbol\xi_T\cdot\boldsymbol n
+\boldsymbol\xi_R\cdot(\boldsymbol r\times\boldsymbol n). \tag{1-3}
\]

因此六维广义法向量为

\[
\boxed{\boldsymbol n_g=
[n_x,n_y,n_z,(\boldsymbol r\times\boldsymbol n)_x,
(\boldsymbol r\times\boldsymbol n)_y,(\boldsymbol r\times\boldsymbol n)_z]}. \tag{1-4}
\]

【源码对照：`Source Code/4.Potential/compute_generalized_normals.m`，第 47–53 行】

```matlab
    nj = zeros(size(centers, 1), 6 * n_bodies); first = 1;
    for b = 1:n_bodies
        idx = first:(first + counts(b) - 1); cols = (b - 1) * 6 + (1:6);
        n_body = normals(idx, :); cg = reshape(body_list{b}.cg, 1, 3);
        nj(idx, cols(1:3)) = n_body;
        nj(idx, cols(4:6)) = cross(centers(idx, :) - cg, n_body, 2);
        first = idx(end) + 1;
    end
```

### 1.3 多体全局排序

若共有 $N_b$ 个浮体，则全局自由度数 $N_d=6N_b$。第 $k$ 个浮体对应的全局索引为

\[
\boxed{\mathcal I_k=6(k-1)+(1:6)}. \tag{1-5}
\]

所有质量、静水恢复、附加质量与阻尼子块均按此顺序嵌入全局矩阵。

【源码对照：`Source Code/5.Force/assemble_mass_matrix.m`，第 26–32 行】

```matlab
    n = numel(mass_props); mass_matrix = zeros(6 * n);
    for b = 1:n
        idx = (b - 1) * 6 + (1:6); m = mass_props(b).mass; inertia = mass_props(b).inertia;
        if ~isequal(size(inertia), [3, 3]) || m <= 0 || any(eig((inertia + inertia.') / 2) <= 0)
            error('CRESTU:MassProperties', 'Body %d mass/inertia are invalid.', b);
        end
        mass_matrix(idx, idx) = blkdiag(m * eye(3), (inertia + inertia.') / 2);
    end
```

### 1.4 复幅值、相位与压力

任一一阶物理量取实部：

\[
f^{(1)}(\boldsymbol x,t)=\Re\{\hat f(\boldsymbol x)e^{\mathrm i\omega t}\}.
\]

线性非定常 Bernoulli 方程为 $p^{(1)}=-\rho\partial_t\phi^{(1)}$，故复压力幅值满足

\[
\boxed{\hat p^{(1)}=-\mathrm i\omega\rho\,\Phi}. \tag{1-6}
\]

两个复幅值的周期平均恒等式为

\[
\overline{\Re\{ae^{\mathrm i\omega t}\}\Re\{be^{\mathrm i\omega t}\}}
=\frac12\Re\{a b^*\}. \tag{1-7}
\]

【源码对照：`Source Code/5.Force/compute_wave_excitation.m`，第 41–42 行】

```matlab
    pressure = -1i * omega * rho * (phi_incident + phi_diffraction);
    force = nj.' * (pressure .* areas);
```

---

## 2. 线性势流控制方程与边界条件

### 2.1 从质量守恒到 Laplace 方程

不可压缩连续性方程为 $\nabla\cdot\boldsymbol u=0$；无旋条件 $\nabla\times\boldsymbol u=0$ 允许定义 $\boldsymbol u=\nabla\phi$。逐项代入可得

\[
0=\nabla\cdot\boldsymbol u
=\nabla\cdot\nabla\phi
=\frac{\partial^2\phi}{\partial x^2}
+\frac{\partial^2\phi}{\partial y^2}
+\frac{\partial^2\phi}{\partial z^2},
\qquad
\boxed{\nabla^2\Phi=0}. \tag{2-1}
\]

源码不在体域内离散 $\nabla^2$，而通过满足 Laplace 方程的 Rankine 基本解将其等价转化为边界积分。

【源码对照：`Source Code/4.Potential/assemble_rankine_matrix.m`，第 74–82 行】

```matlab
            for im = 1:numel(image_weights)
                panel = reshape(image_vertices(j, :, :, im), 4, 3);
                normal = reshape(image_normals(j, :, im), 1, 3);
                [G, dGdn] = rankine_panel_integrals(panel, x, normal);
                Gsum = Gsum + image_weights(im) * G;
                Dsum = Dsum + image_weights(im) * dGdn;
            end
            Sij = Gsum * inv4pi; Dij = -Dsum * inv4pi;
            if i == j, Dij = Dij + 0.5; end
```

### 2.2 势函数分解

总一阶势可写为

\[
\Phi=\Phi_0+\Phi_7+\sum_{j=1}^{6N_b}\xi_j\Phi_j^{(u)}, \tag{2-2}
\]

其中 $\Phi_0$ 为入射势，$\Phi_7$ 为绕射势，$\Phi_j^{(u)}$ 是按单位速度归一化的辐射势。代码的辐射入口按**单位位移**求势，记作 $\Phi_j^{(d)}$，两者关系为

\[
\boxed{\Phi_j^{(d)}=\mathrm i\omega\Phi_j^{(u)}}. \tag{2-3}
\]

【源码对照：`Source Code/4.Potential/solve_radiation_freq.m`，第 38–42 行】

```matlab
    nj = compute_generalized_normals(geom.centers(1:nb, :), geom.normals(1:nb, :), domain.body_list);
    [K, S] = assemble_rankine_matrix(geom.total_panels, geom.centers, geom.normals, geom.vertices, ...
        domain.stats, omega, g, domain.cfg);
    phi = solve_complex_system(K, S * (1i * omega * nj)); phi_body = phi(1:nb, :);
    [A, B, diagnostics] = compute_hydrodynamic_coeffs(phi_body, nj, geom.areas(1:nb), omega, rho);
```

### 2.3 浮体物面边界条件

物面不穿透条件是流体法向速度等于物面法向速度。单位速度归一化时

\[
\boxed{\frac{\partial\Phi_j^{(u)}}{\partial n}=n_j},
\qquad j=1,\ldots,6N_b. \tag{2-4}
\]

【源码对照：`Source Code/4.Potential/compute_generalized_normals.m`，第 49–52 行；代码先构造单位速度归一化所需的 $n_j$】

```matlab
    idx = first:(first + counts(b) - 1); cols = (b - 1) * 6 + (1:6);
    n_body = normals(idx, :); cg = reshape(body_list{b}.cg, 1, 3);
    nj(idx, cols(1:3)) = n_body;
    nj(idx, cols(4:6)) = cross(centers(idx, :) - cg, n_body, 2);
```

按单位位移归一化并采用 $e^{\mathrm i\omega t}$ 后，$\dot\xi_j=\mathrm i\omega\xi_j$，故代码右端为

\[
\boxed{\frac{\partial\Phi_j^{(d)}}{\partial n}=\mathrm i\omega n_j}. \tag{2-5}
\]

【源码对照：`Source Code/4.Potential/solve_radiation_freq.m`，第 39–41 行】

```matlab
    [K, S] = assemble_rankine_matrix(geom.total_panels, geom.centers, geom.normals, geom.vertices, ...
        domain.stats, omega, g, domain.cfg);
    phi = solve_complex_system(K, S * (1i * omega * nj)); phi_body = phi(1:nb, :);
```

固定物体绕射条件理论上为

\[
\boxed{\frac{\partial\Phi_7}{\partial n}=-\frac{\partial\Phi_0}{\partial n}}. \tag{2-6}
\]

当前目录只显式计算右侧入射法向速度，并由载荷函数接收外部给定的 `phi_diffraction`；未找到独立绕射求解入口。

【源码对照：`Source Code/4.Potential/compute_incident_wave.m` 第 47–50 行；`Source Code/5.Force/compute_wave_excitation.m` 第 35–42 行】

```matlab
    phi_I = coefficient * vertical .* phase;
    u = (-1i * k * cos(beta)) * phi_I; v = (-1i * k * sin(beta)) * phi_I;
    w = coefficient * vertical_dz .* phase;
    dphi_I_dn = u .* normals(:, 1) + v .* normals(:, 2) + w .* normals(:, 3);

    if isempty(phi_incident), phi_incident = zeros(size(phi_diffraction), 'like', phi_diffraction); end
    if isempty(phi_diffraction), phi_diffraction = zeros(size(phi_incident), 'like', phi_incident); end
    pressure = -1i * omega * rho * (phi_incident + phi_diffraction);
    force = nj.' * (pressure .* areas);
```

### 2.4 线性自由水面条件与数值阻尼层

在 $z=0$ 处，线性动力学条件与运动学条件分别为

\[
\mathrm i\omega\Phi+g\zeta=0,
\qquad
\mathrm i\omega\zeta=\frac{\partial\Phi}{\partial z}.
\]

由第一式得 $\zeta=-\mathrm i\omega\Phi/g$，代入第二式：

\[
\frac{\partial\Phi}{\partial z}
=\mathrm i\omega\left(-\frac{\mathrm i\omega}{g}\Phi\right)
=\frac{\omega^2}{g}\Phi,
\qquad
\boxed{\partial_z\Phi-\nu\Phi=0,\quad \nu=\omega^2/g}. \tag{2-7}
\]

【源码对照：`Source Code/4.Potential/assemble_rankine_matrix.m`，第 85–86 行；`nu` 是自由面 Robin 系数】

```matlab
    elseif j <= n_body + n_fs
        rowA(j) = Dij - nu(j) * Sij;
```

代码先由色散关系得到 $k_0$，再在 $r>r_{\rm in}$ 区域设置二次阻尼

\[
\mu(r)=\operatorname{clip}\!\left[
\mu_{\rm eff}\left(\frac{r-r_{\rm in}}{r_{\rm out}-r_{\rm in}}\right)^2,
0,\mu_{\rm eff}\right],
\quad
\boxed{\nu_c=k_0(1-\mathrm i\mu)}. \tag{2-8}
\]

【源码对照：`Source Code/4.Potential/assemble_rankine_matrix.m`，第 54–64、85–86 行】

```matlab
    nu = zeros(N, 1); [k0, ~] = solve_wave_dispersion(omega, g, cfg.water_depth);
    sponge_width = cfg.fs.r_outer - cfg.fs.r_inner;
    effective_mu0 = min(2.5, max(cfg.fs.mu0, 12 / (k0 * sponge_width)));
    fs_idx = n_body + (1:n_fs);
    if n_fs > 0
        radii = sqrt(sum(centers(fs_idx, 1:2) .^ 2, 2)); damping = zeros(n_fs, 1);
        mask = radii > cfg.fs.r_inner;
        damping(mask) = effective_mu0 * ((radii(mask) - cfg.fs.r_inner) / sponge_width) .^ 2;
        damping = max(0, min(effective_mu0, damping));
        nu(fs_idx) = k0 * (1 - 1i * damping);
    end
    % ...
    rowA(j) = Dij - nu(j) * Sij;
```

### 2.5 海底与远场条件

有限水深不可穿透海底条件为 $\partial_n\Phi=0$（水平海底等价于 $\partial_z\Phi=0$）；无限水深条件为 $\nabla\Phi\to0$（$z\to-\infty$）。代码对海底仅保留双层项 `Dij`，即齐次 Neumann 条件。

\[
\boxed{\partial_n\Phi=0\quad(z=-h)}. \tag{2-9}
\]

【源码对照：`Source Code/4.Potential/assemble_rankine_matrix.m`，第 85–88 行】

```matlab
            elseif j <= n_body + n_fs
                rowA(j) = Dij - nu(j) * Sij;
            elseif j <= n_body + n_fs + n_sb
                rowA(j) = Dij;
```

在 $e^{\mathrm i\omega t}$ 约定下，外传圆柱波径向相位是 $e^{-\mathrm ikr}$。若远场面法向指向计算域内部，则 Sommerfeld 条件写成

\[
\boxed{\partial_n\Phi=+\mathrm ik\Phi}. \tag{2-10}
\]

【源码对照：`Source Code/4.Potential/assemble_rankine_matrix.m`，第 89–91 行】

```matlab
            else
                % Far-field normals point inward: dphi/dn=+i*k*phi.
                rowA(j) = Dij - 1i * nu(j) * Sij;
            end
```

---

## 3. 入射波与色散关系

有限水深 Airy 波势设为

\[
\Phi_0=\frac{\mathrm igA}{\omega}
\frac{\cosh k(z+h)}{\cosh kh}
\exp[-\mathrm ik(x\cos\beta+y\sin\beta)]. \tag{3-1}
\]

将其代入自由水面条件 $\partial_z\Phi_0=(\omega^2/g)\Phi_0$，在 $z=0$ 约去公共因子：

\[
k\tanh(kh)=\frac{\omega^2}{g},
\qquad
\boxed{\omega^2=gk\tanh(kh)}. \tag{3-2}
\]

【源码对照：`Source Code/4.Potential/solve_wave_dispersion.m`，第 32–34 行】

```matlab
    for iter = 1:50
        kh = k * depth; t = tanh(kh); f = g * k * t - omega^2;
        df = g * t + g * k * depth / (cosh(kh)^2); candidate = max(k - f / df, eps);
```

无限水深 $kh\to\infty$ 时 $\tanh(kh)\to1$，得到 $k=\omega^2/g$。有限水深采用 Newton 迭代：

\[
f(k)=gk\tanh(kh)-\omega^2,
\quad
f'(k)=g\tanh(kh)+gkh\operatorname{sech}^2(kh),
\quad
\boxed{k_{m+1}=k_m-f(k_m)/f'(k_m)}. \tag{3-3}
\]

【源码对照：`Source Code/4.Potential/solve_wave_dispersion.m`，第 29–39 行】

```matlab
    k = omega^2 / g;
    if depth > 0
        k = max(k, omega / sqrt(g * depth));
        for iter = 1:50
            kh = k * depth; t = tanh(kh); f = g * k * t - omega^2;
            df = g * t + g * k * depth / (cosh(kh)^2); candidate = max(k - f / df, eps);
            if abs(candidate - k) <= 1e-12 * max(1, k), k = candidate; break; end
            k = candidate;
        end
    end
    wavelength = 2 * pi / k;
```

对式 (3-1) 求导，得到

\[
\partial_x\Phi_0=-\mathrm ik\cos\beta\,\Phi_0,
\quad
\partial_y\Phi_0=-\mathrm ik\sin\beta\,\Phi_0,
\quad
\partial_z\Phi_0=\frac{\mathrm igA}{\omega}
\frac{k\sinh k(z+h)}{\cosh kh}e^{-\mathrm iks},
\]

从而 $\partial_n\Phi_0=\nabla\Phi_0\cdot\boldsymbol n$。

【源码对照：`Source Code/4.Potential/compute_incident_wave.m`，第 38–50 行】

```matlab
    k = solve_dispersion(omega, g, depth); beta = beta_deg * pi / 180;
    phase = exp(-1i * k * (centers(:, 1) * cos(beta) + centers(:, 2) * sin(beta)));
    coefficient = 1i * g * amplitude / omega;
    if depth > 0
        vertical = cosh(k * (centers(:, 3) + depth)) / cosh(k * depth);
        vertical_dz = k * sinh(k * (centers(:, 3) + depth)) / cosh(k * depth);
    else
        vertical = exp(k * centers(:, 3)); vertical_dz = k * vertical;
    end
    phi_I = coefficient * vertical .* phase;
    u = (-1i * k * cos(beta)) * phi_I; v = (-1i * k * sin(beta)) * phi_I;
    w = coefficient * vertical_dz .* phase;
    dphi_I_dn = u .* normals(:, 1) + v .* normals(:, 2) + w .* normals(:, 3);
```

---

## 4. Green 恒等式、Rankine BEM 与 Hess–Smith 解析面元

### 4.1 Green 第三公式

取未含 $1/(4\pi)$ 的 Rankine 核 $G_0(P,Q)=1/r$，$r=|P-Q|$。对调和函数 $\Phi$ 应用 Green 第二恒等式，并围绕场点挖去半径 $\varepsilon$ 小球；令 $\varepsilon\to0$ 后，小球积分给出立体角系数 $c(P)$，于是

\[
\boxed{c(P)\Phi(P)=\frac1{4\pi}\int_S
\left[G_0(P,Q)\frac{\partial\Phi}{\partial n_Q}
-\Phi(Q)\frac{\partial G_0}{\partial n_Q}\right]\mathrm dS_Q}. \tag{4-1}
\]

光滑边界上 $c=1/2$。定义离散单层、双层系数

\[
S_{ij}=\frac1{4\pi}\int_{S_j}G_0\,\mathrm dS,
\qquad
D_{ij}=-\frac1{4\pi}\int_{S_j}\partial_{n_Q}G_0\,\mathrm dS
+\frac12\delta_{ij}. \tag{4-2}
\]

则直接 BIE 写为 $\sum_jD_{ij}\Phi_j=\sum_jS_{ij}q_j$。

【源码对照：`Source Code/4.Potential/assemble_rankine_matrix.m`，第 67、77–84 行】

```matlab
    A = complex(zeros(N, N)); S_body = complex(zeros(N, n_body)); inv4pi = 1 / (4 * pi);
    % ...
                [G, dGdn] = rankine_panel_integrals(panel, x, normal);
                Gsum = Gsum + image_weights(im) * G;
                Dsum = Dsum + image_weights(im) * dGdn;
            end
            Sij = Gsum * inv4pi; Dij = -Dsum * inv4pi;
            if i == j, Dij = Dij + 0.5; end
            if j <= n_body
                rowA(j) = Dij; rowS(j) = Sij;
```

### 4.2 面元局部坐标变换

令面元中心为 $\boldsymbol x_c$，局部正交基为 $\boldsymbol e_1,\boldsymbol e_2,\boldsymbol e_3$，其中 $\boldsymbol e_3=\boldsymbol n$。变换矩阵按基向量逐行排列：

\[
\mathbf T=
\begin{bmatrix}\boldsymbol e_1^{\mathsf T}\\
\boldsymbol e_2^{\mathsf T}\\
\boldsymbol e_3^{\mathsf T}\end{bmatrix},
\quad
\boldsymbol p_k'=\mathbf T(\boldsymbol p_k-\boldsymbol x_c),
\quad
\boldsymbol x'=\mathbf T(\boldsymbol x-\boldsymbol x_c). \tag{4-3}
\]

【源码对照：`Source Code/4.Potential/rankine_panel_integrals.m`，第 34–56 行】

```matlab
    p_center = mean(P, 1);
    e3 = normal / norm(normal);
    e1 = P(2, :) - P(1, :);
    e1 = e1 - dot(e1, e3) * e3;
    if norm(e1) < 1e-12
        e1 = P(3, :) - P(2, :);
        e1 = e1 - dot(e1, e3) * e3;
    end
    e1 = e1 / norm(e1);
    e2 = cross(e3, e1);
    T = [e1; e2; e3];
    P_loc = (P - p_center) * T';
    X_loc = (X - p_center) * T';
```

### 4.3 边对数项与势积分

局部面元顶点为 $(\xi_k,\eta_k,0)$，场点为 $(x,y,z)$，闭合索引 $k+1$。定义

\[
\Delta\xi_k=\xi_{k+1}-\xi_k,\quad
\Delta\eta_k=\eta_{k+1}-\eta_k,\quad
e_k=\sqrt{\Delta\xi_k^2+\Delta\eta_k^2},
\]

\[
r_k=\sqrt{(x-\xi_k)^2+(y-\eta_k)^2+z^2},
\quad
L_k=\ln\frac{r_k+r_{k+1}+e_k}{r_k+r_{k+1}-e_k}, \tag{4-4}
\]

\[
\cos\theta_k=\Delta\xi_k/e_k,
\quad
\sin\theta_k=\Delta\eta_k/e_k,
\quad
d_k=(\xi_k-x)\sin\theta_k-(\eta_k-y)\cos\theta_k. \tag{4-5}
\]

沿每条直边先作解析积分，再由 Green 定理把面积积分化为边界和，可得

\[
\boxed{I_G=\int_{S_j}\frac{1}{r}\,\mathrm dS
=\sum_k d_kL_k-zI_D},
\qquad
I_D=\int_{S_j}\partial_{n_Q}(1/r)\,\mathrm dS. \tag{4-6}
\]

【源码对照：`Source Code/4.Potential/rankine_panel_integrals.m`，第 62–95 行】

```matlab
    for k = 1:4
        x1 = xi(k);     y1 = eta(k);
        x2 = xi(k + 1);   y2 = eta(k + 1);
        dx = x2 - x1;
        dy = y2 - y1;
        d12 = sqrt(dx^2 + dy^2);
        if d12 < eps_tol, continue; end
        r1 = sqrt((x - x1)^2 + (y - y1)^2 + z^2);
        r2 = sqrt((x - x2)^2 + (y - y2)^2 + z^2);
        cos_th = dx / d12;
        sin_th = dy / d12;
        t0 = (x1 - x) * sin_th - (y1 - y) * cos_th;
        denom_log = r1 + r2 - d12;
        if denom_log > eps_tol
            log_term = log((r1 + r2 + d12) / denom_log);
            G = G + t0 * log_term;
        end
    end
    % ...
    G = G - z * dGdn;
```

### 4.4 立体角、反正切项与法向导数

三角形相对场点的三个向量为 $\boldsymbol a,\boldsymbol b,\boldsymbol c$。Van Oosterom–Strackee 有向立体角公式为

\[
\boxed{\Omega_\triangle=2\operatorname{atan2}\!\left(
\boldsymbol a\cdot(\boldsymbol b\times\boldsymbol c),
abc+(\boldsymbol a\cdot\boldsymbol b)c
+(\boldsymbol b\cdot\boldsymbol c)a
+(\boldsymbol c\cdot\boldsymbol a)b\right)}. \tag{4-7}
\]

【源码对照：`Source Code/4.Potential/rankine_panel_integrals.m`，第 126–129 行】

```matlab
    numerator = dot(a, cross(b, c));
    denominator = norm(a) * norm(b) * norm(c) + dot(a, b) * norm(c) + ...
        dot(b, c) * norm(a) + dot(c, a) * norm(b);
    omega = 2 * atan2(numerator, denominator);
```

四边形由 $(1,2,3)$ 与 $(1,3,4)$ 两三角形组成，$\Omega=\Omega_{123}+\Omega_{134}$。在代码的源法向与顶点方向约定下

\[
\boxed{I_D=-\Omega}. \tag{4-8}
\]

【源码对照：`Source Code/4.Potential/rankine_panel_integrals.m`，第 89–95、126–129 行】

```matlab
    if abs(z) > eps_tol
        a = P(1, :) - X; b = P(2, :) - X; c = P(3, :) - X; d = P(4, :) - X;
        omega = triangle_solid_angle(a, b, c) + triangle_solid_angle(a, c, d);
        dGdn = -omega;
    end
    G = G - z * dGdn;

    numerator = dot(a, cross(b, c));
    denominator = norm(a) * norm(b) * norm(c) + dot(a, b) * norm(c) + ...
        dot(b, c) * norm(a) + dot(c, a) * norm(b);
    omega = 2 * atan2(numerator, denominator);
```

### 4.5 Hess–Smith 单位源面元诱导速度

对 $I_G$ 关于场点坐标求梯度，切向分量化为边对数和，法向分量化为立体角：

\[
u'=\frac1{4\pi}\sum_k\frac{\Delta\eta_k}{e_k}L_k,
\quad
v'=-\frac1{4\pi}\sum_k\frac{\Delta\xi_k}{e_k}L_k,
\quad
w'=-\frac{\Omega}{4\pi}. \tag{4-9}
\]

由于 $[u',v',w']$ 是行向量，返回全局分量为

\[
\boxed{[U,V,W]=[u',v',w']\mathbf T}. \tag{4-10}
\]

【源码对照：`Source Code/3.HessSmith/hess_smith_panel_velocity.m`，第 61–101 行】

```matlab
    for k = 1:4
        x1 = xi(k); y1 = eta(k); x2 = xi(k + 1); y2 = eta(k + 1);
        dx = x2 - x1; dy = y2 - y1; edge_len = hypot(dx, dy);
        if edge_len < eps_tol, continue; end
        r1 = sqrt((x - x1)^2 + (y - y1)^2 + z^2);
        r2 = sqrt((x - x2)^2 + (y - y2)^2 + z^2);
        denom = r1 + r2 - edge_len; numer = r1 + r2 + edge_len;
        if denom > eps_tol && numer > denom
            log_term = log(numer / denom);
            u_loc = u_loc + (dy / edge_len) * log_term;
            v_loc = v_loc - (dx / edge_len) * log_term;
        end
    end
    omega = triangle_solid_angle(q1, q2, q3) + triangle_solid_angle(q1, q3, q4);
    w_loc = -omega;
    vel_glob = ([u_loc, v_loc, w_loc] / (4.0 * pi)) * T;
    u = vel_glob(1); v = vel_glob(2); w = vel_glob(3);
```

### 4.6 经典源强方程与本项目直接 BIE 的区别

经典间接 Hess–Smith 法令

\[
\Phi(P)=\sum_j\frac{\sigma_j}{4\pi}\int_{S_j}\frac1r\,\mathrm dS,
\]

在物面配点取法向导数，得到理论源强方程

\[
\boxed{\sum_j A^{(\sigma)}_{ij}\sigma_j=\mathrm{RHS}_i,
\quad
A^{(\sigma)}_{ij}=\boldsymbol n_i\cdot\boldsymbol v_j(P_i)}. \tag{4-11}
\]

【源码对照：`Source Code/3.HessSmith/hess_smith_panel_velocity.m`，第 98–101 行；该函数提供 $\boldsymbol v_j(P_i)$，但当前主装配未用它求 $\sigma$】

```matlab
    vel_glob = ([u_loc, v_loc, w_loc] / (4.0 * pi)) * T;
    u = vel_glob(1);
    v = vel_glob(2);
    w = vel_glob(3);
```

当前主装配未调用上述速度核，而是组装 $D\Phi=S q$。代入自由面、海底、远场边界条件后，实际未知量方程为

\[
\boxed{\mathbf K\boldsymbol\Phi=\mathbf S_B\boldsymbol q_B}, \tag{4-12}
\]

其中物面列 $K_{ij}=D_{ij}$，自由面列 $D_{ij}-\nu_jS_{ij}$，海底列 $D_{ij}$，远场列 $D_{ij}-\mathrm i kS_{ij}$。

【源码对照：`Source Code/4.Potential/assemble_rankine_matrix.m`，第 70–95 行】

```matlab
    for i = 1:N
        rowA = complex(zeros(1, N)); rowS = complex(zeros(1, n_body)); x = centers(i, :);
        for j = 1:N
            % Gsum and Dsum are accumulated from analytic panel integrals.
            Sij = Gsum * inv4pi; Dij = -Dsum * inv4pi;
            if i == j, Dij = Dij + 0.5; end
            if j <= n_body
                rowA(j) = Dij; rowS(j) = Sij;
            elseif j <= n_body + n_fs
                rowA(j) = Dij - nu(j) * Sij;
            elseif j <= n_body + n_fs + n_sb
                rowA(j) = Dij;
            else
                rowA(j) = Dij - 1i * nu(j) * Sij;
            end
        end
        A(i, :) = rowA; S_body(i, :) = rowS;
    end
```

---

## 5. 一阶水动力系数与载荷

### 5.1 附加质量与辐射阻尼

对单位速度辐射势，定义

\[
\boxed{A_{ij}-\frac{\mathrm i}{\omega}B_{ij}
=\rho\int_{S_0}\Phi_j^{(u)}n_i\,\mathrm dS}. \tag{5-1}
\]

【源码对照：`Source Code/5.Force/compute_hydrodynamic_coeffs.m`，第 42、52 行；代码势为单位位移归一化，后续按式 (5-2) 转换】

```matlab
    potential_integrals = nj.' * (phi_radiation .* areas);
    radiation_force = -1i * omega * rho * potential_integrals;
```

代码求得的是 $\Phi_j^{(d)}=\mathrm i\omega\Phi_j^{(u)}$，因此压力积分为

\[
F_{ij}^{(d)}=-\mathrm i\omega\rho\int\Phi_j^{(d)}n_i\,\mathrm dS
=\omega^2A_{ij}-\mathrm i\omega B_{ij}. \tag{5-2}
\]

逐取实部、虚部即

\[
\boxed{A_{ij}=\Re F_{ij}^{(d)}/\omega^2,
\qquad B_{ij}=-\Im F_{ij}^{(d)}/\omega}. \tag{5-3}
\]

【源码对照：`Source Code/5.Force/compute_hydrodynamic_coeffs.m`，第 42–54 行】

```matlab
    potential_integrals = nj.' * (phi_radiation .* areas);
    symmetry_weights = ones(size(potential_integrals));
    if ~isempty(varargin)
        symmetry = varargin{1};
        if isstruct(symmetry) && isfield(symmetry, 'mode_parity')
            symmetry_weights = symmetry_force_weights(symmetry.mode_parity, ...
                symmetry.mode_parity, symmetry.isx, symmetry.isy);
            potential_integrals = potential_integrals .* symmetry_weights;
        end
    end
    radiation_force = -1i * omega * rho * potential_integrals;
    added_mass_raw = real(radiation_force) / (omega^2);
    damping_raw = -imag(radiation_force) / omega;
```

势流互易关系要求 $A=A^{\mathsf T}$、$B=B^{\mathsf T}$。代码先保留原始非对称残差，再作最近 Frobenius 范数对称投影：

\[
\boxed{A\leftarrow(A+A^{\mathsf T})/2,
\qquad B\leftarrow(B+B^{\mathsf T})/2}. \tag{5-4}
\]

【源码对照：`Source Code/5.Force/compute_hydrodynamic_coeffs.m`，第 55–63 行】

```matlab
    scaleA = max(norm(added_mass_raw, 'fro'), eps); scaleB = max(norm(damping_raw, 'fro'), eps);
    raw_added_mass_symmetry_error = norm(added_mass_raw - added_mass_raw.', 'fro') / scaleA;
    raw_damping_symmetry_error = norm(damping_raw - damping_raw.', 'fro') / scaleB;
    added_mass = 0.5 * (added_mass_raw + added_mass_raw.');
    damping = 0.5 * (damping_raw + damping_raw.');
```

### 5.2 一阶波浪激励力

总固定体散射势为 $\Phi_s=\Phi_0+\Phi_7$。由式 (1-6) 与虚功原理，广义激励力为

\[
\boxed{F_i^{(1)}=-\mathrm i\omega\rho
\int_{S_0}(\Phi_0+\Phi_7)n_i\,\mathrm dS}. \tag{5-5}
\]

其中 $\Phi_0$ 贡献 Froude–Krylov 力，$\Phi_7$ 贡献绕射力。

【源码对照：`Source Code/5.Force/compute_wave_excitation.m`，第 35–42 行】

```matlab
    if isempty(phi_incident), phi_incident = zeros(size(phi_diffraction), 'like', phi_diffraction); end
    if isempty(phi_diffraction), phi_diffraction = zeros(size(phi_incident), 'like', phi_incident); end
    if ~isequal(size(phi_incident), size(phi_diffraction)) || ...
            size(phi_incident, 1) ~= size(nj, 1) || numel(areas) ~= size(nj, 1)
        error('CRESTU:ExcitationShape', 'Incident/diffraction potentials, normals, and areas are inconsistent.');
    end
    pressure = -1i * omega * rho * (phi_incident + phi_diffraction);
    force = nj.' * (pressure .* areas);
```

### 5.3 Haskind 互易校核与能量阻尼

令 $\Psi_j=\Phi_j^{(d)}/(\mathrm i\omega)=\Phi_j^{(u)}$，则 Haskind 型面积积分在代码符号下为

\[
\boxed{F_j^{H}=-\mathrm i\omega\rho m_s
\int_{S_0}\left(n_j\Phi_0+\Psi_j\partial_n\Phi_0\right)\mathrm dS}, \tag{5-6}
\]

其中 $m_s$ 是对称域倍数。

【源码对照：`Source Code/5.Force/compute_haskind_excitation.m`，第 39–46 行】

```matlab
    psi = phi_radiation / (1i * omega);
    for j = 1:ndof
        parity = mode_parity(j, :);
        for h = 1:nh
            [phi_I, dphi_I] = decompose_incident_wave_symmetry(centers, normals, omega, cfg.grav, ...
                cfg.water_depth, headings(h), 1, cfg.isx, cfg.isy, parity);
            integral = sum((nj(:, j) .* phi_I + psi(:, j) .* dphi_I) .* areas);
            force(j, h) = -1i * omega * cfg.rho * cfg.symmetry.multiplicity * integral;
        end
    end
```

角向离散能量关系给出另一套阻尼估计

\[
\boxed{\mathbf B_E=\Re\left\{
\frac{\omega k}{4\pi\rho g^2D_h}
\sum_m\boldsymbol F_H(\theta_m)\boldsymbol F_H(\theta_m)^{\!H}\Delta\theta
\right\}}, \tag{5-7}
\]

有限水深因子 $D_h=\tanh kh+kh\operatorname{sech}^2kh$，无限水深 $D_h=1$。

【源码对照：`Source Code/5.Force/compute_radiation_damping_energy.m`，第 37–49 行】

```matlab
    headings = (0:n_theta - 1) * (360 / n_theta); dtheta = 2 * pi / n_theta;
    force = compute_haskind_excitation(phi_radiation, nj, centers, normals, areas, ...
        omega, cfg, headings, mode_parity);
    [wave_number, ~] = solve_wave_dispersion(omega, cfg.grav, cfg.water_depth);
    if cfg.water_depth > 0
        kh = wave_number * cfg.water_depth;
        depth_factor = tanh(kh) + kh / (cosh(kh)^2);
    else
        depth_factor = 1;
    end
    prefactor = omega * wave_number / (4 * pi * cfg.rho * cfg.grav^2 * depth_factor);
    damping = real(prefactor * (force * force') * dtheta);
    damping = 0.5 * (damping + damping.');
```

### 5.4 静水恢复矩阵

水线面面积、一次矩与二次矩定义为

\[
A_{wp}=\int_{A_{wp}}\mathrm dA,
\quad Q_x=\int x\,\mathrm dA,
\quad Q_y=\int y\,\mathrm dA,
\quad I_{xx}=\int y^2\mathrm dA,
\quad I_{yy}=\int x^2\mathrm dA,
\quad I_{xy}=\int xy\mathrm dA. \tag{5-8}
\]

每个浮体的非零恢复项为

\[
\boxed{
\begin{aligned}
C_{33}&=\rho gA_{wp},&C_{34}&=\rho gQ_y,&C_{35}&=-\rho gQ_x,\\
C_{44}&=\rho gI_{xx}+\rho g\nabla z_B-mgz_G,&
C_{55}&=\rho gI_{yy}+\rho g\nabla z_B-mgz_G,&
C_{45}&=-\rho gI_{xy}.
\end{aligned}} \tag{5-9}
\]

【源码对照：`Source Code/5.Force/compute_hydrostatic_matrix.m`，第 47–56 行】

```matlab
        zb = hs.center_of_buoyancy(3); mass = cfg.mass_props(b).mass; block = zeros(6);
        rg = cfg.rho * cfg.grav;
        block(3, 3) = rg * Awp;
        block(3, 4) = rg * Qy; block(4, 3) = block(3, 4);
        block(3, 5) = -rg * Qx; block(5, 3) = block(3, 5);
        vertical_term = rg * volume * zb - mass * cfg.grav * cg(3);
        block(4, 4) = rg * Ixx + vertical_term;
        block(5, 5) = rg * Iyy + vertical_term;
        block(4, 5) = -rg * Ixy; block(5, 4) = block(4, 5);
        idx = (b - 1) * 6 + (1:6); C(idx, idx) = block;
```

### 5.5 频域刚体方程与 RAO

时域线性方程

\[
(\mathbf M+\mathbf A)\ddot{\boldsymbol\xi}
+\mathbf B\dot{\boldsymbol\xi}+\mathbf C\boldsymbol\xi
=\boldsymbol F_{exc}
\]

代入 $\boldsymbol\xi(t)=\Re\{\hat{\boldsymbol\xi}e^{\mathrm i\omega t}\}$，使用 $\dot\xi=\mathrm i\omega\xi$、$\ddot\xi=-\omega^2\xi$，逐项得到

\[
\boxed{\mathbf Z(\omega)\hat{\boldsymbol\xi}
=\hat{\boldsymbol F}_{exc},\quad
\mathbf Z=-\omega^2(\mathbf M+\mathbf A)+\mathrm i\omega\mathbf B+\mathbf C}, \tag{5-10}
\]

【源码对照：`Source Code/5.Force/solve_rao.m`，第 44–45 行】

```matlab
    w = omegas(k); Z = -w^2 * (M + added_mass(:, :, k)) + 1i * w * damping(:, :, k) + hydrostatic;
    dynamic_matrix(:, :, k) = Z; rconds(k) = rcond(Z);
```

故运动响应传递函数（单位波幅 RAO）为

\[
\boxed{\hat{\boldsymbol\xi}=\mathbf Z^{-1}\hat{\boldsymbol F}_{exc}}. \tag{5-11}
\]

【源码对照：`Source Code/5.Force/solve_rao.m`，第 41–50 行】

```matlab
    nh = size(excitation, 2); M = assemble_mass_matrix(cfg.mass_props);
    xi = complex(zeros(ndof, nh, nf)); dynamic_matrix = complex(zeros(ndof, ndof, nf)); rconds = zeros(nf, 1);
    for k = 1:nf
        w = omegas(k); Z = -w^2 * (M + added_mass(:, :, k)) + 1i * w * damping(:, :, k) + hydrostatic;
        dynamic_matrix(:, :, k) = Z; rconds(k) = rcond(Z);
        if rconds(k) < 1e-12, warning('CRESTU:IllConditionedRAO', 'RAO matrix at omega=%g has rcond=%g.', w, rconds(k)); end
        xi(:, :, k) = Z \ excitation(:, :, k);
    end
    rao = struct('omega', omegas, 'xi', xi, 'amplitude', abs(xi), 'phase_deg', angle(xi) * 180 / pi, ...
        'mass_matrix', M, 'dynamic_matrix', dynamic_matrix, 'rcond', rconds, ...
        'convention', 'exp(+i*omega*t)');
```

---

## 6. Pinkster 近场平均漂移力

### 6.1 二阶周期平均与总分解

对同频复幅值，二次积的平均含 $1/2$。完整近场平均载荷写成

\[
\boxed{\overline{\boldsymbol F}^{(2)}
=\boldsymbol F_{d1}+\boldsymbol F_{d2}+\boldsymbol F_{d3}+\boldsymbol F_{d4}}. \tag{6-1}
\]

源码字段映射为：`term_waterline` $\leftrightarrow F_{d1}$；`term_quadratic_velocity` $\leftrightarrow F_{d2}$；`term_translation_gradient` $\leftrightarrow$ 裁剪的 $F_{d3}$；`term_rotation_force` $\leftrightarrow$ 裁剪的 $F_{d4}$。源码数组相加顺序不影响总和。

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 91–94、130–139 行】

```matlab
    term_waterline = zeros(6, heading_count);
    term_quadratic_velocity = zeros(6, heading_count);
    term_rotation_force = zeros(6, heading_count);
    term_translation_gradient = zeros(6, heading_count);
    % ...
    total = term_waterline + term_quadratic_velocity + term_rotation_force + term_translation_gradient;
    drift = struct( ...
        'term_waterline', term_waterline, ...
        'term_quadratic_velocity', term_quadratic_velocity, ...
        'term_rotation_force', term_rotation_force, ...
        'term_translation_gradient', term_translation_gradient, ...
        'total', total, ...
```

### 6.2 第一项：相对波面升高的水线积分

线性自由面升高由动力学条件得

\[
\hat\zeta=-\frac{\mathrm i\omega}{g}\Phi_{WL}. \tag{6-2}
\]

水线点 $\boldsymbol r=(x-x_G,y-y_G,0)$ 的刚体竖向位移为

\[
\hat z_B=\xi_3+(\boldsymbol\xi_R\times\boldsymbol r)_z
=\xi_3+(y-y_G)\xi_4-(x-x_G)\xi_5. \tag{6-3}
\]

因此相对波高

\[
\boxed{\hat\zeta_r=\hat\zeta-
[\xi_3+(y-y_G)\xi_4-(x-x_G)\xi_5]}. \tag{6-4}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 100–104 行】

```matlab
    wave_elevation = -1i * state.omega * waterline_potential(:, heading_index) / cfg.grav;
    body_elevation = translation(3) ...
        + rotation(1) * (nodes(:, 2) - center_of_gravity(2)) ...
        - rotation(2) * (nodes(:, 1) - center_of_gravity(1));
    relative_elevation = wave_elevation - body_elevation;
```

周期平均 $\overline{\zeta_r^2}=|\hat\zeta_r|^2/2$，水静压力三角形积分给出

\[
\boxed{\boldsymbol F_{d1}=\frac12\rho g
\oint_{WL}|\hat\zeta_r|^2\boldsymbol n_g^{WL}\,\mathrm dl}. \tag{6-5}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 74–78、100–107 行】

```matlab
    midpoint = 0.5 * (nodes + nodes(next_node, :));
    line_normal = [segment(:, 2), -segment(:, 1)] ./ max(segment_length, eps);
    line_position = [midpoint, zeros(node_count, 1)] - center_of_gravity;
    line_normal_three = [line_normal, zeros(node_count, 1)];
    line_generalized_normal = [line_normal_three, cross(line_position, line_normal_three, 2)];

    wave_elevation = -1i * state.omega * waterline_potential(:, heading_index) / cfg.grav;
    body_elevation = translation(3) ...
        + rotation(1) * (nodes(:, 2) - center_of_gravity(2)) ...
        - rotation(2) * (nodes(:, 1) - center_of_gravity(1));
    relative_elevation = wave_elevation - body_elevation;
    midpoint_elevation = 0.5 * (relative_elevation + relative_elevation(next_node));
    term_waterline(:, heading_index) = 0.5 * cfg.rho * cfg.grav ...
             * (line_generalized_normal.' * (abs(midpoint_elevation) .^ 2 .* segment_length));
```

变量映射：$\Phi_{WL}\leftrightarrow$ `waterline_potential`，$\zeta\leftrightarrow$ `wave_elevation`，$z_B\leftrightarrow$ `body_elevation`，$\zeta_r\leftrightarrow$ `relative_elevation`，$\mathrm dl\leftrightarrow$ `segment_length`。

### 6.3 第二项：Bernoulli 速度平方项

非定常 Bernoulli 方程保留到二阶：

\[
p=-\rho\left(\partial_t\phi+gz+\frac12|\nabla\phi|^2\right).
\]

一阶速度复幅值为 $\hat{\boldsymbol u}=\nabla\Phi$，且

\[
\overline{|\boldsymbol u^{(1)}|^2}
=\frac12\hat{\boldsymbol u}\cdot\hat{\boldsymbol u}^{*}.
\]

于是源码采用的广义二阶动压项为

\[
\boxed{\boldsymbol F_{d2}=-\frac12\rho
\int_{S_0}(\nabla\Phi\cdot\nabla\Phi^*)\boldsymbol n_g\,\mathrm dS}. \tag{6-6}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 109–112 行】

```matlab
    surface_velocity = reshape(velocity(:, :, heading_index), panel_count, 3);
    speed_squared = real(sum(surface_velocity .* conj(surface_velocity), 2));
    term_quadratic_velocity(:, heading_index) = -0.5 * cfg.rho ...
             * (generalized_normal.' * (speed_squared .* areas));
```

变量映射：$\nabla\Phi\leftrightarrow$ `surface_velocity`，$|\nabla\Phi|^2\leftrightarrow$ `speed_squared`，$\boldsymbol n_g\leftrightarrow$ `generalized_normal`。

### 6.4 第三项：位移—压力梯度耦合

把一阶压力从平均物面 $\boldsymbol x$ Taylor 展开到位移后的物面 $\boldsymbol x+\boldsymbol X$：

\[
p^{(1)}(\boldsymbol x+\boldsymbol X)
=p^{(1)}(\boldsymbol x)+\boldsymbol X\cdot\nabla p^{(1)}(\boldsymbol x)+O(\varepsilon^3).
\]

故完整理论项为

\[
\boxed{\boldsymbol F_{d3}^{\rm full}
=-\frac12\Re\int_{S_0}
[\boldsymbol X^*\cdot\nabla\hat p^{(1)}]\boldsymbol n_g\,\mathrm dS,
\quad \nabla\hat p^{(1)}=-\mathrm i\omega\rho\nabla\Phi}. \tag{6-7}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 121–124 行；这是当前最接近的 Hessian 收缩切片，但缺少完整式中的显式 $-\mathrm i\omega$ 因子与角位移】

```matlab
    potential_hessian = reshape(hessian(panel_index, :, :, heading_index), 3, 3);
    directional_gradient = potential_hessian * conj(translation(:));
    translated_normal_acceleration(panel_index) = real( ...
        dot(directional_gradient, normals(panel_index, :).'));
```

当前代码仅取平移 $\boldsymbol\xi_T$，用势 Hessian $\mathbf H_\Phi=\nabla\nabla\Phi$ 构造

\[
\boxed{\boldsymbol F_{d3}^{\rm code}
=-\frac12\rho\int_{S_0}
\Re\{\boldsymbol n\cdot\mathbf H_\Phi\boldsymbol\xi_T^*\}
\boldsymbol n_g\,\mathrm dS}. \tag{6-8}
\]

式 (6-8) 与式 (6-7) 并不代数等价：源码中没有显式 $-\mathrm i\omega$ 压力梯度因子，也未把角位移纳入 $\boldsymbol X$。因此此项必须解释为“平移—势 Hessian 裁剪项”，不能宣称是完整 Pinkster 压力梯度项。

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 118–127 行】

```matlab
    translated_normal_acceleration = zeros(panel_count, 1);
    for panel_index = 1:panel_count
        potential_hessian = reshape(hessian(panel_index, :, :, heading_index), 3, 3);
        directional_gradient = potential_hessian * conj(translation(:));
        translated_normal_acceleration(panel_index) = real( ...
            dot(directional_gradient, normals(panel_index, :).'));
    end
    term_translation_gradient(:, heading_index) = -0.5 * cfg.rho ...
             * (generalized_normal.' * (translated_normal_acceleration .* areas));
```

变量映射：$\mathbf H_\Phi\leftrightarrow$ `potential_hessian`，$\boldsymbol\xi_T^*\leftrightarrow$ `conj(translation)`，法向收缩量 $\leftrightarrow$ `translated_normal_acceleration`。

### 6.5 第四项：角运动旋转耦合

一阶力随小转角 $\boldsymbol\Omega$ 旋转到惯性系时，向量增量为 $\boldsymbol\Omega\times\boldsymbol F^{(1)}$。周期平均完整结构可写为

\[
\boxed{\boldsymbol F_{d4}^{\rm full}
=\frac12\Re\{\boldsymbol\Omega^*\times\boldsymbol F^{(1)}\}
+\boldsymbol F_{\rm inertia/moment}^{(2)}}. \tag{6-9}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 114–116 行；仅第一项的平移力部分存在，惯性/力矩项未实现】

```matlab
    % The requested cross term uses the first-order excitation force only.
    term_rotation_force(1:3, heading_index) = 0.5 * real(cross( ...
        conj(rotation(:)), state.first_order_force(1:3, heading_index)));
```

当前实现只写入前三个平移分量：

\[
\boxed{\boldsymbol F_{d4}^{\rm code}(1{:}3)
=\frac12\Re\{\boldsymbol\Omega^*\times\boldsymbol F_T^{(1)}\}},
\qquad
\boldsymbol F_{d4}^{\rm code}(4{:}6)=\boldsymbol0. \tag{6-10}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_nearfield.m`，第 114–116 行】

```matlab
    % The requested cross term uses the first-order excitation force only.
    term_rotation_force(1:3, heading_index) = 0.5 * real(cross( ...
        conj(rotation(:)), state.first_order_force(1:3, heading_index)));
```

变量映射：$\boldsymbol\Omega\leftrightarrow$ `rotation`，$\boldsymbol F_T^{(1)}\leftrightarrow$ `state.first_order_force(1:3,...)`。惯性力矩耦合未见源码对应变量或公式。

### 6.6 表面速度与 Hessian 的数值重构

在每个面元局部切平面，以邻点坐标 $(u_m,v_m)$ 拟合二次多项式

\[
\Delta\Phi_m\approx
a_1u_m+a_2v_m+\frac12a_3u_m^2+a_4u_mv_m+\frac12a_5v_m^2. \tag{6-11}
\]

采用高斯权 $w_m=\exp[-4(r_m/R)^2]$ 与 Tikhonov 正则化，法方程为

\[
\boxed{\boldsymbol a=(\mathbf A^H\mathbf W\mathbf A+\lambda s\mathbf I)^{-1}
\mathbf A^H\mathbf W\Delta\boldsymbol\Phi}. \tag{6-12}
\]

【源码对照：`Source Code/6.MeanDriftLoads/estimate_surface_kinematics.m`，第 65–69、89–99 行】

```matlab
    radius = sqrt(local_u .^ 2 + local_v .^ 2);
    support_radius = max(radius(end), eps);
    weight = exp(-4 * (radius / support_radius) .^ 2);
    design = [local_u, local_v, 0.5 * local_u .^ 2, local_u .* local_v, 0.5 * local_v .^ 2];
    weighted_design = design .* sqrt(weight);
    % ...
    potential_difference = phi(neighbor_index, heading_index) - phi(panel_index, heading_index);
    weighted_rhs = potential_difference .* sqrt(weight);
    normal_matrix = weighted_design' * weighted_design;
    scale = max(trace(real(normal_matrix)) / 5, 1);
    coefficient = (normal_matrix + regularization * scale * eye(5)) ...
            \ (weighted_design' * weighted_rhs);
    velocity(panel_index, :, heading_index) = ...
        coefficient(1) * tangent_one(panel_index, :) ...
        + coefficient(2) * tangent_two(panel_index, :) ...
        + dphi_dn(panel_index, heading_index) * normals(panel_index, :);
```

### 6.7 无量纲平均漂移系数

平移力采用参考尺度 $\tfrac12\rho gA^2L$，力矩再乘 $L$：

\[
\boxed{C_{d,i}=\frac{F_i^{(2)}}{\tfrac12\rho gA^2L},\ i=1,2,3;
\qquad
C_{d,i}=\frac{M_i^{(2)}}{\tfrac12\rho gA^2L^2},\ i=4,5,6}. \tag{6-13}
\]

【源码对照：`Source Code/6.MeanDriftLoads/export_drift_loads.m`，第 34–35、52–54 行】

```matlab
    force_scale = 0.5 * rho * g * wave_amplitude^2 * L;
    moment_scale = force_scale * L;
    % ...
    scale = force_scale * ones(rows, 1); scale(dof > 3) = moment_scale;
    output = table(frequency, heading, method, body, dof, load_value, load_value ./ scale, ...
        'VariableNames', {'omega_rad_s', 'heading_deg', 'method', 'body', 'dof', 'load', 'Cd'});
```

---

## 7. Maruo–Newman 远场漂移力（源码扩展模块）

对角向试验函数 $\psi_\theta$，代码用 Green 恒等式构造 Kochin 函数

\[
\boxed{H(\theta)=-\int_{S_0}
[\psi_\theta\partial_n\Phi-\Phi\partial_n\psi_\theta]\,\mathrm dS}. \tag{7-1}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_farfield.m`，第 43–56 行】

```matlab
    phase = exp(-1i * k * (centers(:, 1) * ct + centers(:, 2) * st));
    if cfg.water_depth > 0
        fz = cosh(k * (centers(:, 3) + cfg.water_depth)) / cosh(k * cfg.water_depth);
        dfz = k * sinh(k * (centers(:, 3) + cfg.water_depth)) / cosh(k * cfg.water_depth);
    else
        fz = exp(k * centers(:, 3)); dfz = k * fz;
    end
    psi = fz .* phase;
    dpsi = psi .* (-1i * k * (normals(:, 1) * ct + normals(:, 2) * st)) + dfz .* phase .* normals(:, 3);
    H(t, h) = -sum((psi .* state.dphi_dn(:, h) - state.phi(:, h) .* dpsi) .* areas);
```

角导数采用中心差分 $H'_m=(H_{m+1}-H_{m-1})/(2\Delta\theta)$，再按动量通量求 surge、sway、yaw。源码仅输出整体系统的 1、2、6 分量，不分配到各浮体。

\[
\boxed{H'(\theta_m)\approx\frac{H_{m+1}-H_{m-1}}{2\Delta\theta}}. \tag{7-2}
\]

【源码对照：`Source Code/6.MeanDriftLoads/compute_drift_farfield.m`，第 59–76 行】

```matlab
    Hp = (circshift(H, -1, 1) - circshift(H, 1, 1)) / (2 * dtheta);
    loads = zeros(6, nh); denominator = nu + cfg.water_depth * (k^2 - nu^2);
    if cfg.water_depth <= 0, denominator = nu; end
    prefactor = cfg.rho * k^3 / (8 * pi * nu) * (nu / denominator);
    cross_prefactor = cfg.rho * k * state.omega / (2 * nu);
    loads(1, h) = real(prefactor * sum(cos(theta) .* H(:, h) .* conj(H(:, h))) * dtheta + ...
        cross_prefactor * cos(beta) * interp_periodic(theta, H(:, h), pi + beta));
    loads(2, h) = real(prefactor * sum(sin(theta) .* H(:, h) .* conj(H(:, h))) * dtheta + ...
        cross_prefactor * sin(beta) * interp_periodic(theta, H(:, h), pi + beta));
    yaw_integral = sum(Hp(:, h) .* conj(H(:, h))) * dtheta;
    loads(6, h) = real(-1i * cfg.rho * k^2 / (8 * pi * nu) * (nu / denominator) * yaw_integral + ...
        1i * cfg.rho * state.omega / (4 * nu) * yaw_incident_pair);
```

---

## 8. 网格几何、面积、排水体积与浮心

四边形面元以两条对角线 $\boldsymbol d_{13}$、$\boldsymbol d_{24}$ 近似面积向量：

\[
\boxed{S_p=\frac12|\boldsymbol d_{13}\times\boldsymbol d_{24}|,
\qquad \boldsymbol n_p=\frac{\boldsymbol d_{13}\times\boldsymbol d_{24}}
{|\boldsymbol d_{13}\times\boldsymbol d_{24}|}}. \tag{8-1}
\]

【源码对照：`Source Code/2.Mesh/read_bmf.m`，第 78–92 行】

```matlab
        d13 = P(3, :) - P(1, :);
        d24 = P(4, :) - P(2, :);
        n_raw = cross(d13, d24);
        mag = norm(n_raw);
        if mag < 1e-12
            n_raw = cross(P(2, :) - P(1, :), P(4, :) - P(1, :));
            mag = norm(n_raw);
        end
        normals(k, :) = n_raw / mag;
        areas(k) = 0.5 * mag;
```

由散度定理，$\nabla\cdot(x\boldsymbol e_x)=1$，所以闭合体积可由 $\int_Sxn_x\mathrm dS$ 求得；$y,z$ 同理。离散后

\[
\boxed{V_x=\sum_px_pn_{x,p}S_p,
\quad V_y=\sum_py_pn_{y,p}S_p,
\quad V_z=\sum_pz_pn_{z,p}S_p,
\quad V=(V_x+V_y+V_z)/3}. \tag{8-2}
\]

【源码对照：`Source Code/2.Mesh/read_bmf.m`，第 120–128 行】

```matlab
        sym_factor = (1 + isx) * (1 + isy);
        hydro.Vx = sum(centers(:, 1) .* normals(:, 1) .* areas) * sym_factor;
        hydro.Vy = sum(centers(:, 2) .* normals(:, 2) .* areas) * sym_factor;
        hydro.Vz = sum(centers(:, 3) .* normals(:, 3) .* areas) * sym_factor;
        hydro.V_mean = (hydro.Vx + hydro.Vy + hydro.Vz) / 3.0;
        r_dot_n = sum(centers .* normals, 2);
        dV = (1.0 / 3.0) * r_dot_n .* areas;
```

---

## 9. 对称性恢复

若输出自由度在 $x=0,y=0$ 平面的奇偶性为 $p_i^x,p_i^y\in\{\pm1\}$，解势的奇偶性为 $p_j^x,p_j^y$，则缩减域积分恢复权重为

\[
\boxed{w_{ij}=(1+p_i^xp_j^x)^{I_x}(1+p_i^yp_j^y)^{I_y}}, \tag{9-1}
\]

其中未启用的对称面指数取零。奇偶性不相容时权重自动为零。

【源码对照：`Source Code/4.Potential/symmetry_force_weights.m`，第 29–31 行】

```matlab
    weights = ones(size(output_parity, 1), size(solution_parity, 1));
    if isx, weights = weights .* (1 + output_parity(:, 1) * solution_parity(:, 1).'); end
    if isy, weights = weights .* (1 + output_parity(:, 2) * solution_parity(:, 2).'); end
```

入射波的镜像波向及奇偶权重按群投影平均；双对称面时

\[
\boxed{\Phi^{(p_x,p_y)}=\frac14[
\Phi(\beta)+p_x\Phi(180^\circ-\beta)+p_y\Phi(-\beta)
+p_xp_y\Phi(180^\circ+\beta)]}. \tag{9-2}
\]

【源码对照：`Source Code/4.Potential/decompose_incident_wave_symmetry.m`，第 44–56 行】

```matlab
    if isx && isy
        flags = [0, 0;1, 0;0, 1;1, 1];
        headings = [beta_deg, 180 - beta_deg, -beta_deg, 180 + beta_deg];
        weights = [1, parity(1), parity(2), parity(1) * parity(2)];
    end
    phi_component = complex(zeros(size(centers, 1), 1));
    dphi_component = complex(zeros(size(centers, 1), 1));
    for q = 1:size(flags, 1)
        [phi, dphi] = compute_incident_wave(centers, normals, omega, g, depth, headings(q), amplitude);
        phi_component = phi_component + weights(q) * phi;
        dphi_component = dphi_component + weights(q) * dphi;
    end
    scale = 2^(isx + isy); phi_component = phi_component / scale; dphi_component = dphi_component / scale;
```

---

## 10. 源码—理论公式全景对照表

| 公式编号 | 物理意义 | 理论数学表达式 | 对应 `.m` 源码文件及行号/子函数 | 核心英文变量名 |
|---|---|---|---|---|
| (1-1) | 局部—全局刚体变换 | $\boldsymbol x=\boldsymbol x_0+R_z\boldsymbol x_k$ | `2.Mesh/transform_body_mesh.m:29-36` | `Rz`, `x0`, `vertices`, `centers` |
| (1-4) | 六维广义法向 | $[\boldsymbol n,\boldsymbol r\times\boldsymbol n]$ | `4.Potential/compute_generalized_normals.m:47-53` | `nj`, `n_body`, `centers`, `cg` |
| (1-5) | 多体自由度映射 | $6(k-1)+(1:6)$ | `5.Force/assemble_mass_matrix.m:26-32` | `idx`, `mass_matrix` |
| (1-6) | 线性动水压力 | $p=-i\omega\rho\Phi$ | `5.Force/compute_wave_excitation.m:41-42` | `pressure`, `omega`, `rho` |
| (2-1) | Laplace 方程边界化 | $\nabla^2\Phi=0$ | `4.Potential/assemble_rankine_matrix.m:74-82` | `G`, `dGdn`, `Sij`, `Dij` |
| (2-5) | 单位位移辐射条件 | $\partial_n\Phi=i\omega n_j$ | `4.Potential/solve_radiation_freq.m:38-42` | `S`, `omega`, `nj`, `phi` |
| (2-6) | 固定体绕射条件 | $\partial_n\Phi_7=-\partial_n\Phi_0$ | 入射法向量：`4.Potential/compute_incident_wave.m:47-50`；独立绕射求解入口缺失 | `dphi_I_dn`, `phi_diffraction` |
| (2-8) | 自由面阻尼层 | $\nu_c=k(1-i\mu)$ | `4.Potential/assemble_rankine_matrix.m:54-64,85-86` | `nu`, `damping`, `effective_mu0` |
| (2-9) | 海底 Neumann 条件 | $\partial_n\Phi=0$ | `4.Potential/assemble_rankine_matrix.m:87-88` | `rowA`, `Dij` |
| (2-10) | 远场辐射条件 | $\partial_n\Phi=ik\Phi$ | `4.Potential/assemble_rankine_matrix.m:89-91` | `rowA`, `nu`, `Sij` |
| (3-2) | 有限水深色散 | $\omega^2=gk\tanh kh$ | `4.Potential/solve_wave_dispersion.m:29-39` | `k`, `kh`, `t`, `f`, `df` |
| (3-1) | Airy 入射势 | $igA\cosh k(z+h)e^{-iks}/(\omega\cosh kh)$ | `4.Potential/compute_incident_wave.m:38-50` | `phase`, `vertical`, `coefficient`, `phi_I` |
| (4-1) | Green 第三公式 | $c\Phi=(4\pi)^{-1}\int(Gq-\Phi G_n)dS$ | `4.Potential/assemble_rankine_matrix.m:67,77-84` | `Gsum`, `Dsum`, `Sij`, `Dij` |
| (4-3) | 面元局部坐标 | $x'=T(x-x_c)$ | `4.Potential/rankine_panel_integrals.m:34-56` | `e1`, `e2`, `e3`, `T`, `P_loc`, `X_loc` |
| (4-6) | Rankine 势解析积分 | $I_G=\sum d_kL_k-zI_D$ | `4.Potential/rankine_panel_integrals.m:62-95` | `d12`, `r1`, `r2`, `t0`, `log_term`, `G` |
| (4-7) | 三角形有向立体角 | $2\operatorname{atan2}(N,D)$ | `4.Potential/rankine_panel_integrals.m:126-129` | `numerator`, `denominator`, `omega` |
| (4-9) | Hess–Smith 诱导速度 | 边对数切向项与立体角法向项 | `3.HessSmith/hess_smith_panel_velocity.m:61-101` | `u_loc`, `v_loc`, `w_loc`, `vel_glob` |
| (4-11) | 经典间接源强方程 | $A^{(\sigma)}\sigma=RHS$ | 速度核存在：`3.HessSmith/hess_smith_panel_velocity.m`；主装配未调用 | `u`, `v`, `w` |
| (4-12) | 项目直接 BIE | $K\Phi=S_Bq_B$ | `4.Potential/assemble_rankine_matrix.m:70-95` | `A`, `S_body`, `rowA`, `rowS` |
| (5-3) | 附加质量与阻尼 | $A=\Re F/\omega^2,\ B=-\Im F/\omega$ | `5.Force/compute_hydrodynamic_coeffs.m:42-54` | `potential_integrals`, `radiation_force`, `added_mass_raw`, `damping_raw` |
| (5-5) | 一阶激励力 | $-i\omega\rho\int(\Phi_0+\Phi_7)n_i dS$ | `5.Force/compute_wave_excitation.m:35-42` | `phi_incident`, `phi_diffraction`, `pressure`, `force` |
| (5-6) | Haskind 关系 | $-i\omega\rho\int(n_j\Phi_0+\Psi_j\Phi_{0,n})dS$ | `5.Force/compute_haskind_excitation.m:39-46` | `psi`, `phi_I`, `dphi_I`, `integral`, `force` |
| (5-7) | 能量法阻尼 | $B_E\propto\int F_HF_H^H d\theta$ | `5.Force/compute_radiation_damping_energy.m:37-49` | `depth_factor`, `prefactor`, `damping` |
| (5-9) | 静水恢复 | 水线面积/矩与浮力—重力项 | `5.Force/compute_hydrostatic_matrix.m:47-56` | `Awp`, `Qx`, `Qy`, `Ixx`, `Iyy`, `Ixy`, `block` |
| (5-10) | 动力刚度 | $Z=-\omega^2(M+A)+i\omega B+C$ | `5.Force/solve_rao.m:41-47` | `M`, `Z`, `dynamic_matrix`, `xi` |
| (6-5) | Pinkster 水线项 | $\tfrac12\rho g\oint|\zeta_r|^2n_gdl$ | `6.MeanDriftLoads/compute_drift_nearfield.m:74-78,100-107` | `relative_elevation`, `line_generalized_normal`, `segment_length` |
| (6-6) | 速度平方项 | $-\tfrac12\rho\int|\nabla\Phi|^2n_gdS$ | `6.MeanDriftLoads/compute_drift_nearfield.m:109-112` | `surface_velocity`, `speed_squared`, `term_quadratic_velocity` |
| (6-8) | 平移—Hessian 裁剪项 | $-\tfrac12\rho\int\Re[n\cdot H_\Phi\xi_T^*]n_gdS$ | `6.MeanDriftLoads/compute_drift_nearfield.m:118-127` | `potential_hessian`, `directional_gradient`, `term_translation_gradient` |
| (6-10) | 转动—一阶力裁剪项 | $\tfrac12\Re(\Omega^*\times F_T^{(1)})$ | `6.MeanDriftLoads/compute_drift_nearfield.m:114-116` | `rotation`, `first_order_force`, `term_rotation_force` |
| (6-12) | 表面导数重构 | 正则加权最小二乘 | `6.MeanDriftLoads/estimate_surface_kinematics.m:65-69,89-99` | `design`, `weight`, `coefficient`, `velocity` |
| (6-13) | 漂移系数 | $F/(0.5\rho gA^2L)$ | `6.MeanDriftLoads/export_drift_loads.m:34-35,52-54` | `force_scale`, `moment_scale`, `Cd` |
| (7-1) | Kochin 函数 | $-\int(\psi\Phi_n-\Phi\psi_n)dS$ | `6.MeanDriftLoads/compute_drift_farfield.m:43-56` | `psi`, `dpsi`, `H` |
| (8-1) | 面元面积与法向 | $S=|d_{13}\times d_{24}|/2$ | `2.Mesh/read_bmf.m:78-92` | `d13`, `d24`, `n_raw`, `areas` |
| (8-2) | 排水体积 | $V_x=\sum xn_xS$ 等 | `2.Mesh/read_bmf.m:120-128` | `Vx`, `Vy`, `Vz`, `V_mean`, `dV` |
| (9-1) | 对称积分权重 | $(1+p_i^xp_j^x)(1+p_i^yp_j^y)$ | `4.Potential/symmetry_force_weights.m:29-31` | `output_parity`, `solution_parity`, `weights` |

---

## 11. 实现一致性核查清单

1. **相位：** 入射相位使用 `exp(-1i*k*s)`，压力使用 `-1i*omega*rho*phi`，RAO 动力刚度使用 `+1i*omega*B`，三者与 $e^{+i\omega t}$ 一致。
2. **法向：** 主 BIE 的物面法向按网格存储方向使用；远场代码明确假设法向指向计算域内部。改变顶点顺序会同时改变 $D_{ij}$、广义法向和载荷符号。
3. **归一化：** 辐射求解入口是单位位移归一化；公式 (5-1) 是单位速度归一化。二者必须通过 $\Phi^{(d)}=i\omega\Phi^{(u)}$ 转换。
4. **直接法/间接法：** `assemble_rankine_matrix` 的 `A` 是直接 BIE 的混合边界矩阵，不是经典源强影响矩阵 $A^{(\sigma)}$。
5. **绕射闭合：** 当前目录只能消费 `phi_diffraction`，不能单独证明其由式 (2-6) 求得；集成应用必须提供该势或补充求解入口。
6. **近场二阶项：** `term_translation_gradient` 与 `term_rotation_force` 应按式 (6-8)、(6-10) 命名和验证，不应直接等同于完整式 (6-7)、(6-9)。
7. **多体近场：** `compute_drift_nearfield` 返回 6 行载荷，调用时应确认 `mesh` 是单体网格或明确整体广义法向定义；远场实现只给整体 surge、sway、yaw。
8. **数值验证：** 应同时检查 $A,B$ 互易残差、能量法阻尼半正定性、近远场漂移力一致性、网格加密收敛性和海绵层半径敏感性。

---

## 12. 经典参考文献

1. Hess, J. L., and Smith, A. M. O. (1964). Calculation of nonlifting potential flow about arbitrary three-dimensional bodies.
2. Newman, J. N. (1977). *Marine Hydrodynamics*. MIT Press.
3. Faltinsen, O. M. (1990). *Sea Loads on Ships and Offshore Structures*. Cambridge University Press.
4. Pinkster, J. A. (1980). *Low Frequency Second Order Wave Exciting Forces on Floating Structures*. Delft University of Technology.
5. Newman, J. N. (1974). Second-order, slowly-varying forces on vessels in irregular waves.
