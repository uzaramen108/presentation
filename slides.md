---
theme: seriph
title: 유체-막 상호작용 시뮬레이션
info: |
  ## NS + Brinkman + Darcy + ALE
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
---

<div class="absolute inset-0 bg-white z-0" />
<img src="./images/logo.png"
  class="absolute inset-0 w-full h-full object-contain opacity-15 z-0 pointer-events-none" />

<div class="relative z-10 flex flex-col items-center justify-center h-full" color='black'>

# 유체-막 상호작용 시뮬레이션

<div class="text-xl opacity-80 mb-8 mt-4">
  Navier-Stokes + Brinkman + Darcy + ALE
</div>

</div>

<div class="absolute bottom-6 right-6 text-xl z-10" color='black'>
  윤현준 / 박기성 / 박찬서 / 이민용
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 1. Develop하려는 문제 상황(응용1)

<div class="text-xl opacity-80 mb-6">반투막이 있는 원형 관로에서의 유체 거동</div>

<div grid="~ cols-2 gap-8" class="text-left">

<div v-click>

### 🔹 물리적 설정

원형 단면의 3D 파이프 내부를 유체가 흐르다가 중간에 위치한
반투막(semi-permeable membrane)과 만나는 문제입니다.

<br>

<div class="-mt-4 text-xs opacity-80 bg-[rgb(210,220,220)] rounded-md px-2 py-1 leading-tight">

<style scoped>
table { 
  margin: 0 !important;
  width: 100%;
}
th, td { 
  padding: 1px 6px !important; /* 상하 패딩을 1px로 극한까지 줄임 */
  line-height: 1.1 !important;
  border-bottom-width: 1px !important;
}
</style>

| 파라미터 | 값 | 단위 |
|---|---|---|
| 파이프 길이 $L$ | 5.0 | m |
| 파이프 반지름 $R$ | 0.5 | m |
| 막 위치 $x_m$ | 3.75 | m |
| 유체 밀도 $\rho$ | 1.0 | kg/m³ |
| 유체 점성 $\mu$ | 0.01 | Pa·s |
| 막 투과율 $\kappa$ | 5×10⁻³ | m² |
| 막 장력 $T$ | 25 | N/m |

</div>

</div>

<div v-click>

### 🔹 입구 조건

Poiseuille 포물선 프로파일에 사인파 시간 변화를 적용합니다:

$$u_{inlet}(r,t) = u_{max} \sin\!\left(\frac{\pi t}{2T}\right)\left(1 - \frac{r^2}{R^2}\right)$$

### 🔹 계산 목표

<div class="mt-2 text-sm opacity-80">

- **속도장** $u(x,t)$: 막 전후 유체의 흐름 패턴
- **압력장** $p(x,t)$: 막 전후 압력 분포 및 점프
- **막 변위** $w(x,t)$: 압력 하중에 의한 막의 휨

</div>

<div class="w-full flex justify-center mt-4">
  <img src="./images/screenshot1.png"
    class="w-full max-h-[140px] object-contain rounded shadow-sm border border-gray-200/50" />
</div>

</div>
</div>
---
layout: center
class: text-center
---

# The PDE Problem: Navier-Stokes
<div class="text-xl opacity-80 mb-8">
  비압축성 유체의 지배 방정식
</div>

<div grid="~ cols-2 gap-8" class="text-left mt-8">

<div>

### 1. Momentum & Continuity
유체의 운동량 보존과 비압축성을 나타냅니다.

$$
\rho \left( \frac{\partial u}{\partial t} + u \cdot \nabla u \right) = \nabla \cdot \sigma(u, p) + f
$$
$$
\nabla \cdot u = 0
$$
- $u$: 속도 (Velocity)
- $p$: 압력 (Pressure)
- $f$: 외력 (Body force)

</div>

<div v-click>

### 2. Stress Tensor
뉴턴 유체의 응력 텐서를 정의합니다.

$$
\sigma(u, p) = 2\mu\epsilon(u) - pI
$$
$$
\epsilon(u) = \frac{1}{2}(\nabla u + (\nabla u)^T)
$$
- $2\mu\epsilon$: 유체의 변형에 의한 점성 마찰력
- $-pI$: 압력에 의한 수직 응력

</div>
</div>

---
layout: center
---

# Step 1: 임시 속도 $u^*$ 구하기 (1/3)
<div class="text-xl opacity-80 mb-8">기본 지배 방정식과 시간 이산화</div>

<div class="text-sm mt-8">

### 1. 원래의 운동량 보존 방정식 (강형, Strong Form)
가장 기본이 되는 나비에-스토크스 방정식입니다.
$$ \rho \frac{\partial u}{\partial t} + \rho (u \cdot \nabla u) = \nabla \cdot (2\mu\epsilon(u)) - \nabla p + f $$

<br>

### 2. 시간 이산화 (Time Discretization)
컴퓨터가 계산할 수 있도록 시간 미분 $\frac{\partial u}{\partial t}$를 $\frac{u^* - u^n}{\Delta t}$ (코드의 `(u - u_n) / k`)로 쪼갭니다.
이때, 점성항에는 현재와 과거의 평균을 쓰는 **Crank-Nicolson (CN) 기법**을 적용하여 수치적 안정성을 높입니다.
$$ \rho \frac{u^* - u^n}{\Delta t} + \rho (u \cdot \nabla u) - \nabla \cdot \left(2\mu\epsilon\left(\frac{u^* + u^n}{2}\right)\right) + \nabla p^n = f $$

</div>

---
layout: center
---

# Step 1: 임시 속도 $u^*$ 구하기 (2/3)
<div class="text-xl opacity-80 mb-8">대류항 선형화 및 부분 적분(약형 변환)</div>

<div class="text-sm mt-8">

### 3. 대류항의 선형화 (Linearization of Convection)
비선형항인 $u \cdot \nabla u$를 해결해야 합니다. 속도 벡터 두 개가 곱해져 있어 방정식이 안 풀리므로, 앞의 $u$는 과거 데이터를 이용해 추정치(외삽)를 만들고, 뒤의 $\nabla u$는 미분항에 CN 기법을 씁니다.
* **앞부분 (Adams-Bashforth 외삽)**: $u_{AB} = \frac{3}{2}u^n - \frac{1}{2}u^{n-1}$
* **뒷부분 (Crank-Nicolson 미분)**: $\nabla u_{CN} = \nabla\left(\frac{u^* + u^n}{2}\right)$

<br>

### 4. 부분 적분 (Integration by Parts)
양변에 테스트 함수 $v$를 곱하고 전체 영역에 대해 적분합니다. 이때 2계 미분인 점성항과 1계 미분인 압력항의 미분 차수를 낮춰줍니다.
* **점성항 변환**: $\int -\nabla \cdot (2\mu\epsilon) \cdot v \Rightarrow \int \mu \nabla u : \nabla v$
* **압력항 변환**: $\int \nabla p \cdot v \Rightarrow \int -p (\nabla \cdot v)$

</div>

---
layout: two-cols
layoutClass: gap-8
---

# Step 1: 임시 속도 $u^*$ 구하기 (3/3)
<div class="text-xl opacity-80 mb-8">최종 수식과 코드 매칭</div>

::left::

<div class="text-[13px] mt-4">

### 5. 최종 수식 ($F_1 = 0$)
앞선 시간 이산화, 선형화, 약형 변환을 모두 합치면 코드와 완벽히 일치하는 최종 잔차(Residual) 방정식이 됩니다.

$$
\begin{aligned}
F_1 &= \int \rho \frac{u^* - u^n}{\Delta t} \cdot v
&+ \int \rho (u_{AB} \cdot \nabla u_{CN}) \cdot v
&+ \int \mu \nabla u_{CN} : \nabla v
\end{aligned}
$$
$$
\begin{aligned}
&- \int p^n (\nabla \cdot v)
&+ \int f \cdot v = 0
\end{aligned}
$$

(여기서 $u^*$는 미지수 $u$를 의미합니다.)

</div>

::right::

<div class="overflow-y-auto max-h-[450px] shadow-lg rounded-md border border-gray-200/20 text-sm mt-12">

```python {all} twoslash
# Step 1: Tentative velocity step (u_s)

# 1. 시간 미분항
F1 = rho / k * dot(u - u_n, v) * dx

# 2. 대류항 (Adams-Bashforth + Crank-Nicolson)
F1 += inner(dot(1.5 * u_n - 0.5 * u_n1, 
                0.5 * nabla_grad(u + u_n)), v) * dx

# 3. 점성항 (부분적분 완료, CN 적용)
# 4. 압력항 (부분적분 완료, 이전 스텝 p_ 사용)
F1 += 0.5 * mu * inner(grad(u + u_n), grad(v)) * dx \
      - dot(p_, div(v)) * dx

# 5. 외력항 (f = 0)
F1 += dot(f, v) * dx

a1 = form(lhs(F1))
L1 = form(rhs(F1))
```
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 2-3. 약형 도출


<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div class="text-[13px] mt-4" v-click>

### 1. 수식 변환 과정
진짜 속도와 가짜 속도의 차이는 "새롭게 업데이트될 압력의 구배(기울기)" 때문에 발생합니다.
$$ \rho \frac{u^{n+1} - u^*}{\Delta t} = -\nabla (p^{n+1} - p^n) $$

양변에 발산($\nabla \cdot$)을 취하고, $\nabla \cdot u^{n+1} = 0$을 대입하면 **poisson equation** 형태가 됩니다. (여기서 $\phi = p^{n+1} - p^n$ 이라 둡니다).
$$ \nabla^2 \phi = \frac{\rho}{\Delta t} \nabla \cdot u^* $$

### 2. 약형 (Weak Form) 변환
테스트 함수 $q$를 곱하고 좌변에 부분 적분을 적용합니다.
$$ \int \nabla \phi \cdot \nabla q = -\int \frac{\rho}{\Delta t} (\nabla \cdot u^*) q $$

</div>

<div class="text-[13px] mt-4" v-click>

### 1. 속도 보정 수식 (Velocity Projection)
Step 2에서 구한 식을 다시 $u^{n+1}$에 대해 정리합니다.
$$ u^{n+1} = u^* - \frac{\Delta t}{\rho} \nabla \phi $$

### 2. 약형 (Weak Form) 변환
전체 식에 밀도 $\rho$와 테스트 함수 $v$를 곱하여 적분형(질량 행렬 형태)으로 만듭니다.
$$ \int \rho u^{n+1} \cdot v = \int \rho u^* \cdot v - \int \Delta t \nabla \phi \cdot v $$

- 좌변: 구하려는 진짜 속도 $u^{n+1}$ (코드의 `u`)
- 우변 첫 번째: 이미 구한 가짜 속도 $u^*$ (코드의 `u_s`)
- 우변 두 번째: 이미 구한 압력 증분의 구배 $\nabla \phi$

이 단계를 통과한 $u^{n+1}$은 완벽하게 비압축성 조건을 만족하게 됩니다.

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 2-3. 약형 적용 및 차수 할당

<br>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 Step 3: 속도 보정

$$\int_\Omega \rho\, u^{n+1}\cdot v = \int_\Omega \rho\, u^*\cdot v - \int_\Omega \Delta t\,\nabla\phi\cdot v$$

이 단계를 거친 $u^{n+1}$은 **$\nabla\cdot u = 0$을 정확히 만족**합니다.

$$p^{n+1} = p^n + \phi$$

```python
# NS Step 3: velocity correction
a3_ufl = form(rho_v * dot(u_t, vf) * dxf)
L3_ufl = form(rho_v * dot(u_, vf) * dxf
            - dt * dot(nabla_grad(p_ - p_n), vf) * dxf)
```

</div>

<div v-click>

### 🔹 함수 공간: Taylor-Hood P2/P1

속도-압력 쌍은 **inf-sup 조건(LBB)** 을 만족해야 발진하지 않습니다.

-비압축성 유체 계산에서 속도($u$)와 압력($p$)에 똑같은 차수의 함수(예: 둘 다 P1)를 사용하면, 계산 결과에서 압력이 마치 체스판 무늬처럼 위아래로 진동(Oscillation)이 일어남.

<div style="background-color:rgb(220,220,210)">

| 변수 | 차수 | 이유 |
|---|---|---|
| 속도 $u$ | P2 | inf-sup 조건(압력보다 큰 차수) |
| 압력 $p$ | P1 |  |

</div>
</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 3-1. 메시 설정: 파이프 단독 구조

<div class="text-xl opacity-80 mb-6">Brinkman 방식 — 분리된 2개 메시</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 메시 구성 (초기)
```python
# 유체 메시: 3D 원기둥
cylinder = gmsh.model.occ.addCylinder(
    0, 0, 0, L, 0, 0, R)

# 막 메시: 별도 2D 디스크
disk = gmsh.model.occ.addDisk(0, 0, 0, R, R)
gmsh.model.occ.rotate(...)
gmsh.model.occ.translate(..., mem_pos, 0, 0)
```

두 메시가 **완전히 분리**되어 있습니다.

</div>

<div v-click>

### 🔹 분리된 메시의 문제점
```
fluid_mesh ──┐
              ├── 서로 다른 메시
mem_mesh   ──┘

```

→ 내부 경계면 공유 불가<br>
→ dS 표면 적분 불가<br>
→ FacetNormal을 dx 안에 사용 시 오류<br>
→ KDTree 매핑 필요
</div>

</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 3-1. 메시 설정: 파이프 단독 구조

<div class="text-xl opacity-80 mb-6">Brinkman 방식 — 분리된 2개 메시</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 경계 조건

**유체:**
- 입구: Poiseuille 속도 프로파일
- 출구: 압력 $p = 0$
- 벽면: No-slip ($u = 0$)

**막:**
- 테두리: 고정단 ($w = 0$)

</div>

<div v-click>


### 🔹 압력-막 연결 (KDTree 매핑)
```python
# 유체 압력 DOF 좌표
q_coords = Q_f.tabulate_dof_coordinates()
# 막 DOF 좌표
w_coords = W_m.tabulate_dof_coordinates()

# 최근접 매핑
tree = cKDTree(q_coords)
_, W_to_Qf = tree.query(w_coords)
```

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 3-3. 메시 설정: Fragment 통합 구조

<div class="text-xl opacity-80 mb-6">Darcy + ALE 방식 — 내부 경계면으로 통합</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 Gmsh Fragment 연산

파이프를 막 위치에서 **두 조각으로 나누고** 막 디스크가 공유 경계면이 되도록 합니다:
```python
    c1 = gmsh.model.occ.addCylinder(0,   0, 0, x_m,     0, 0, R)
    c2 = gmsh.model.occ.addCylinder(x_m, 0, 0, L - x_m, 0, 0, R)
    dk = gmsh.model.occ.addDisk(x_m, 0, 0, R, R)
    gmsh.model.occ.rotate([(2, dk)], x_m, 0, 0, 0, 1, 0, np.pi / 2)
    out_tags, out_map = gmsh.model.occ.fragment(
        [(3, c1), (3, c2)], [(2, dk)]
    )
```
<div class="mt-4 text-xs opacity-80 bg-[rgb(200,200,240)] rounded-md px-2 py-1 leading-tight">

<style scoped>
table { 
  margin: 0 !important;
  width: 100%;
}
th, td { 
  padding: 6px 6px !important; /* 상하 패딩을 1px로 극한까지 줄임 */
  line-height: 1.4 !important;
  border-bottom-width: 2px !important;
}
</style>

| Physical Group | fluid (vol) | inlet | outlet | membrane |
| :--- | :--- | :--- | :--- | :--- |
| **역할** | 3D 유체 도메인 | 입구 경계 | 출구 경계 | **내부 경계면** |
</div>
</div>

<div v-click>

### 🔹 Fragment 후의 장점
```
기존:                     Fragment 후:
두 개 메시               하나의 통합 메시
KDTree 매핑              geom_map으로 정확한 인덱스
dS 불가                  dS_m (tag=4) 직접 적분 가능
FacetNormal 오류         FacetNormal 사용 가능
```

### 🔹 2D Submesh 추출

막 계산을 위한 2D submesh를 3D 메시에서 추출합니다:
```python
sub_m, emap_m, vmap_m, geom_map_raw = create_submesh(
    msh, msh.topology.dim - 1, ft.find(4)
)
geom_map = np.array(geom_map_raw, dtype=np.int32)

# 동기화: sub_m.geometry.x[i] = msh.geometry.x[geom_map[i]]
# → 매 스텝 msh가 이동하면 sub_m도 함께 갱신
```

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 4. Brinkman 저항 모델

<div class="text-xl opacity-80 mb-6">막을 다공성 매질로 근사</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 Brinkman 방정식

NS에 **부피 저항항** $\frac{\mu}{\kappa}u$를 추가한 형태입니다:

$$\rho\frac{\partial u}{\partial t} + \rho(u\cdot\nabla u) = -\nabla p + \mu\nabla^2 u - \underbrace{\frac{\mu}{\kappa}u}_{\text{Brinkman 저항}}$$

- $\kappa$: 막의 투과율 (작을수록 저항 증가)
- $\mu/\kappa$: 물리적 저항 계수

### 🔹 NS 방정식의 위계
```
Navier-Stokes (완전)
      ↓ Re → 0
  Stokes
      ↓ 다공성 매질 균질화 + 점성항 유지
  Brinkman         ← 반투막에 적용
      ↓ 점성항 추가 무시
  Darcy 법칙
```

</div>

<div v-click>

### 🔹 초기 구현 vs 개선된 구현

**① 초기 (문제 있음):**

$$R_m = 50\text{ (임의)}, \quad |x - x_m| < 0.2\text{m (두께 0.4m)}$$
```
α
50──────────────
   ████████████
───────────────
   ←  0.4m  →  너무 두꺼움
```

**② 개선 (tanh 집중화):**

$$\alpha(x) = \frac{\mu}{\kappa} \cdot \frac{1}{2}\!\left(1 + \tanh\frac{t_{half} - |x - x_m|}{\varepsilon}\right)$$
```
α
│    █
μ/κ  ███
│   █████
───────────
   ←0.04m→  얇고 집중
```

$t_{half} \approx 0.02\text{m}$, $\alpha_{max} = \mu/\kappa$ (물리량 기반)

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 4-1. Brinkman가 작용하는 방식

<div class="text-xl opacity-80 mb-6">압력 분포 및 물리적 타당성</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 기대되는 물리 현상

* 반투막 통과 시 Darcy 법칙에 의해 **압력 강하(불연속)** 가 반드시 발생해야 합니다:

$$\Delta p = p^- - p^+ = \frac{\mu}{\kappa}(u \cdot n) \cdot d_{mem}$$

<br>

```
실제(Darcy):       p⁻ ━━━┤불연속┝━━━ p⁺
Brinkman:   p⁻ ━━━╲_____╱━━━ p⁺  (완만)
```




</div>

<div v-click>


### 🔹 남은 한계

<br>

* Brinkman 항은 **부피 적분(Volume integral)**:

$$\int_{\Omega} \frac{\mu}{\kappa}\cdot\mathbf{1}_{\Gamma_m^\delta}\cdot u\cdot v\;dx \leftarrow \text{3D 영역에 분산}$$

실제 막의 압력 점프는 **표면 적분(Surface integral)** 이어야 합니다:

$$\int_{\Gamma_m} [\![p]\!]\cdot v\cdot n\;dS \leftarrow \text{2D 경계면에 집중}$$

* Brinkman은 x, y, z **모든 방향 동등하게 저항** 을 받으며 막의 모든 방향에서 통과, 하지만 실제 막은 법선 방향으로만 통과

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 5. Darcy 표면항 — 물리적으로 올바른 해법

<div class="text-xl opacity-80 mb-6">막을 2D 내부 경계면으로 처리</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 Darcy 점프 조건

막을 **2D 내부 경계면** $\Gamma_m$으로 정의하면, 막 통과 유동을 Darcy 법칙으로 정확히 기술할 수 있습니다:

$$[\![p]\!] = \frac{\mu}{\kappa}(u\cdot n), \qquad [\![p]\!] = p^- - p^+$$

이를 NS 운동량 약형에 **표면 적분항**으로 추가합니다:

$$F_1 += \int_{\Gamma_m} \frac{\mu}{\kappa} \langle u \rangle \cdot n^{+} \, \langle v \rangle \cdot n^{+} \, dS$$

- $\langle u \rangle = \text{avg}(u)$: 막 양쪽 속도 평균 (수치 안정성)
- $n^+$: 막의 법선 벡터 — **dS에서만 합법, dx에서는 오류**
- $\mu/\kappa$: 투과율에서 유도된 물리량

</div>

<div v-click>

### 🔹 실제 거동 vs 세 가지 모델 비교

<br>

<div style="background-color:rgb(210,220,240)">

| | 실제 거동 | Brinkman | Darcy (표면항) |
|---|:---:|:---:|:---:|
| 압력 점프 | ✅ 불연속 | △ 완만 | ✅ 강제 |
| 저항 계수 | 기준 | △ 근사 | ✅ $\kappa$ |
| 법선 방향만 | ✅ | ❌ 전방향 | ✅ |
| 접선 차단 | ✅ | ❌ | ✅ BJS |
| 적분 형태 | — | 부피 $dx$ | 표면 $dS$ |

</div>
</div>
</div>
---
theme: seriph
class: text-center
highlighter: shiki
---

# 6. ALE 기법의 필요성

<div class="text-xl opacity-80 mb-6">Arbitrary Lagrangian-Eulerian</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 세 관점의 차이

| 관점 | 격자 움직임 | 적합 대상 | 한계 |
|---|---|---|---|
| **Eulerian** | 고정 | 유체 | 움직이는 경계 추적 불가 |
| **Lagrangian** | 물질과 함께 | 고체 | 유체 적용 시 격자 뒤틀림 |
| **ALE** | 독립적 | FSI | — |

### 🔹 왜 ALE가 필요한가?

막이 압력을 받아 변형되면 유체 도메인 $\Omega_f(t)$가 시간에 따라 변합니다.
```
t=0:  ──────[막]──────   (평평)
t=Δt: ────[막]────        (휘어짐)

고정 메시: 변형된 형상을 인식 못함
           → t=0 행렬로 t=T 계산 → 물리 위반 → 발산
```

</div>

<div v-click>

### 🔹 ALE의 역할 분담
```
구조물 경계면(막):
  격자가 막과 함께 이동 (Lagrangian)
            ↓
  유체 내부:
  격자가 꼬이지 않도록 부드럽게 재배치 (Eulerian)
            ↓
  그 사이:
  독립적인 속도 w_mesh로 이동 (ALE)
```

### 🔹 행렬을 매 스텝 재조립해야 하는 이유
```python
# 루프 안에서 매 스텝:
A1.zeroEntries()
assemble_matrix(A1, a1_ufl, bcs=bcu)
A1.assemble()
s1.setOperators(A1)  # 전처리기도 갱신
```

메시 좌표 변경 → $dx$, $FacetNormal$, 야코비안 모두 변경
→ 재조립 없으면 물리 법칙 위반

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 7. ALE 수식: 대류항 수정

<div class="text-xl opacity-80 mb-6">움직이는 격자 위의 Navier-Stokes</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 ALE 대류 속도

격자 자체가 $w_{mesh}$의 속도로 움직이므로, 유체가 느끼는 **상대적 대류 속도**가 달라집니다.

**기존 대류항:**
$$\rho(u\cdot\nabla u)$$

**ALE 대류항 (메시 속도 반영):**
$$\rho\underbrace{(u - w_{mesh})}_{\text{상대 속도}}\cdot\nabla u$$

이 수정은 Donea et al. (2004)의 ALE 방정식 14번 식에 해당합니다.

</div>

<div v-click>

### 🔹 전체 ALE-NS + Darcy 약형

$$\underbrace{\int_\Omega \rho\frac{u^*-u^n}{\Delta t}\cdot v}_{\text{시간 미분}}
+ \underbrace{\int_\Omega \rho(u^n - w_{mesh})\cdot\nabla u^n\cdot v}_{\text{ALE 대류항}}$$

$$+ \underbrace{\int_\Omega \sigma:\epsilon(v)\;dx}_{\text{점성 + 압력}}
+ \underbrace{\int_{\Gamma_m}\frac{\mu}{\kappa}\langle u\rangle \cdot n\;\langle v\rangle \cdot n\;dS}_{\text{Darcy 표면 저항}} = 0$$

### 🔹 메시 속도 계산

$$w_{mesh} = \frac{d_{ALE}}{\Delta t}$$

$d_{ALE}$: 조화 확장으로 구한 메시 변위 (다음 슬라이드)

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 8. 조화 확장 (Harmonic Extension)

<div class="text-xl opacity-80 mb-6">격자 뒤틀림 없이 막 변위를 도메인으로 전파</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 라플라스 방정식

막이 변위 $\Delta w$만큼 움직이면, 내부 격자들이 겹치지 않도록
**라플라스 방정식**으로 변위를 부드럽게 전파합니다:

$$-\nabla\cdot(k\nabla d) = 0 \quad \text{in } \Omega_f$$

경계 조건:

$$d = \Delta w\,\hat{e}_x \;\text{ on } \Gamma_m, \qquad d = 0 \;\text{ on } \partial\Omega_f$$

- $d$: 메시 변위 벡터
- $\Delta w = w^{n+1} - w^n$: 막의 **증분** 변위
- $\hat{e}_x$: x방향 단위벡터 (막이 x방향으로 휨)
- $\delta$: 0 나누기 방지 오프셋

라플라스 방정식의 해는 경계 조건을 만족하면서 **가장 매끄러운 함수**임이 수학적으로 보장됩니다.

</div>

<div v-click>

### 🔹 가변 강성 (Variable Stiffness)

막에 가까울수록 강성 $k$를 크게 설정하여 좁은 격자의 뒤틀림을 방지합니다:

$$k = \frac{1}{(|x - x_m| + \delta)^3}$$
```
파이프 x방향 강성 분포:

입구                막(3.75)               출구
0  ──[유연]──────[뻣뻣]──────[유연]──── 5
                    ↑
              k 최대 (격자 거의 안 움직임)
```

| 위치 | $k$ | 격자 거동 |
|---|---|---|
| 막 바로 옆 | ~400 | 거의 안 움직임 |
| 막에서 0.5m | ~3.3 | 약간 이동 |
| 막에서 2m | ~0.2 | 크게 이동 가능 |

<div class="mt-2 p-2 bg-red-900 bg-opacity-20 rounded text-xs">
⚠️ 격자가 뒤틀리면 야코비안이 음수 → 행렬 특이 → 계산 붕괴 (직접 경험)
</div>

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 9-1. 막 변위 계산 (Poisson)

<div class="text-xl opacity-80 mb-6">압력 점프로부터 막의 변형 계산</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 막 변위 지배 방정식

탄성막의 변위 $w$는 장력 $T$와 압력 하중 $[\![p]\!]$에 의한 **포아송 방정식**으로 기술됩니다:

$$-T\Delta_\Gamma w = [\![p]\!] = p^- - p^+$$

경계 조건: $w = 0$ on $\partial\Gamma_m$ (막 테두리 고정단)

약형:

$$T\int_{\Gamma_m} \nabla w \cdot \nabla v \, dA = \int_{\Gamma_m} [\![p]\!] \cdot v \, dA$$

### 🔹 압력 점프 $[\![p]\!]$ 계산 방법

막 변위를 반영한 위치에서 압력을 샘플링합니다:

$$[\![p]\!]_i = p(x_m + w_i - \varepsilon\mathbf{n}_i) - p(x_m + w_i + \varepsilon\mathbf{n}_i)$$

</div>

<div v-click>

### 🔹 동적 샘플링의 중요성

막이 변형되면 법선 벡터 $\mathbf{n}$이 더 이상 x축과 평행하지 않습니다.
```
초기 (평평한 막):    변형 후 (휜 막):
법선 = (1,0,0)       법선 = 사선 방향
x축 샘플링 = 맞음    x축만 샘플링 = 오차!
```

따라서 **현재 막 법선 방향**으로 샘플링해야 합니다.

$$\varepsilon_{dyn} = \max(\varepsilon_0,\; 1.5 \cdot w_{max} + 10^{-4})$$

막 변위가 클수록 샘플링 오프셋도 동적으로 확장됩니다.

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 9-2. 막 법선 벡터 계산

<div class="text-xl opacity-80 mb-6">비대칭 변형에 대한 일반화된 법선 벡터</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 법선 벡터 유도

막 표면: $\mathbf{F}(y,z) = (x_m + w(y,z),\; y,\; z)$

접선 벡터:
$$T_y = \left(\frac{\partial w}{\partial y}, 1, 0\right), \quad T_z = \left(\frac{\partial w}{\partial z}, 0, 1\right)$$

법선 벡터 ($T_y \times T_z$):
$$\mathbf{n} = \frac{1}{\sqrt{1 + \left(\frac{\partial w}{\partial y}\right)^2 + \left(\frac{\partial w}{\partial z}\right)^2}}\left(1,\; -\frac{\partial w}{\partial y},\; -\frac{\partial w}{\partial z}\right)$$

</div>

<div v-click>

### 🔹 편미분을 최소자승법으로 계산

각 DOF에서 이웃 DOF들을 이용해 $\partial w/\partial y$, $\partial w/\partial z$를 **독립적으로** 추정합니다:

$$\text{minimize} \sum_j \left|\Delta w_j - \frac{\partial w}{\partial y}\Delta y_j - \frac{\partial w}{\partial z}\Delta z_j\right|^2$$

2×2 정규 방정식:

$$\underbrace{\begin{pmatrix} \sum \Delta y^2 & \sum \Delta y\Delta z \\ \sum \Delta y\Delta z & \sum \Delta z^2 \end{pmatrix}}_{A^TA} \begin{pmatrix} \partial w/\partial y \\ \partial w/\partial z \end{pmatrix} = \underbrace{\begin{pmatrix} \sum \Delta w\Delta y \\ \sum \Delta w\Delta z \end{pmatrix}}_{A^Tb}$$

### 🔹 r 대칭 가정과의 차이

| | r 대칭 가정 (기존) | 최소자승법 (수정) |
|---|---|---|
| 비대칭 변형 | ❌ 오차 발생 | ✅ 정확 |
| r=0 특이점 | ❌ 0 나누기 | ✅ 없음 |
| 비대칭 장애물 | ❌ 불가 | ✅ 가능 |

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 10. 전체 코드 흐름

<div class="text-xl opacity-80 mb-6">매 타임스텝의 연성 계산 순서</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 스태거드 결합 순서
```
타임스텝 n → n+1:

① 입구 속도 갱신 u_in.interpolate(iv)
       ↓
② NS Step 1: u* 계산
   (행렬 재조립 + Darcy dS_m 포함)
       ↓
③ NS Step 2: 압력 p^{n+1}
       ↓
④ NS Step 3: 속도 u^{n+1}
       ↓
⑤ 막 양단 압력 점프 계산
   pressure_jump_dynamic()
   (법선 방향 + MPI Allreduce)
       ↓
⑥ 막 변위 Poisson: w^{n+1}
       ↓
⑦ ALE 조화 확장: d_ale
   (가변 강성, 증분 경계 조건)
       ↓
⑧ 메시 좌표 갱신 + sub_m 동기화
   msh.geometry.x += d_inc
   sub_m.geometry.x = msh.geometry.x[geom_map]
       ↓
⑨ bb_tree 갱신 후 반복
```

</div>

<div v-click>

### 🔹 주요 함수 구조
```python
# 막 법선 방향 압력 점프 (최소자승법)
def pressure_jump_dynamic(p_fn, w_arr,
        mem_dof_coords, gtree, eps):
    # ① ∂w/∂y, ∂w/∂z 최소자승 계산
    # ② n = (1, -∂w/∂y, -∂w/∂z) / |n|
    # ③ 법선 방향 ±ε 샘플링
    # ④ MPI Allreduce (이중 계산 방지)
    return p_minus - p_plus

# KSP 솔브 헬퍼
def ksp_solve(solver, b_vec, x_fn,
              L_f, a_f=None, bcs=None):
    b_vec.localForm().set(0)
    assemble_vector(b_vec, L_f)
    if a_f and bcs:
        apply_lifting(b_vec, [a_f], [bcs])
    b_vec.ghostUpdate(...)
    if bcs: set_bc(b_vec, bcs)
    solver.solve(b_vec, x_fn.x.petsc_vec)
    x_fn.x.scatter_forward()
```

</div>
</div>

---
layout: two-cols
layoutClass: gap-8
---

<script setup>
import { ref } from 'vue'
const isExpanded = ref(false)
</script>

# 10-1. 전체 코드 보기
ALE + Darcy 통합 코드

::left::

<div class="text-sm mt-4">

### 🔹 핵심 파라미터
<div class="text-xs opacity-80 mt-2">

| 항목 | 값 |
|---|---|
| 메시 크기 (최대) | 0.12 m |
| 메시 크기 (최소) | 0.025 m |
| 타임스텝 수 | 2,000 |
| 시뮬레이션 시간 | 2.0 s |
| $\Delta t$ | 1×10⁻³ s |

</div>

### 🔹 출력 파일
- `results/u.bp`: 속도장
- `results/p.bp`: 압력장
- `results/w.bp`: 막 변위
- `results/d.bp`: ALE 변위

</div>

::right::

<button @click="isExpanded = true" class="mt-8 px-4 py-2 bg-gray-800 text-white text-sm rounded shadow-md hover:bg-gray-700 transition-all flex items-center gap-2">
  <carbon:zoom-in /> 코드 전체화면으로 보기
</button>

<div :class="isExpanded ? 'fixed inset-4 z-50 bg-white shadow-2xl rounded-xl p-8 overflow-y-auto border border-gray-300' : 'overflow-y-auto max-h-[450px] shadow-lg rounded-md border border-gray-200/20 text-sm mt-4'">

  <div v-show="isExpanded" class="flex justify-between items-center mb-4 border-b pb-4">
    <div class="text-xl font-bold text-gray-800">FSI ALE + Darcy 전체 코드</div>
    <button @click="isExpanded = false" class="px-4 py-1.5 bg-red-500 text-white rounded hover:bg-red-600 transition-all text-sm flex items-center gap-1">
      <carbon:close /> 닫기
    </button>
  </div>
```python {all} twoslash
  # 코드를 여기에 붙여넣으세요
```

</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 11. 시뮬레이션 결과 분석

<div class="text-xl opacity-80 mb-6">속도장 · 압력장 · 막 변위 비교</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 속도장 (Velocity Field)

<div class="w-full flex flex-col items-center mt-2">
  <img src="./images/image1.png"
    class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">
    막 위치($x=3.75$)에서 속도 감소 확인
  </div>
</div>

</div>

<div v-click>

### 🔹 압력장 비교

<div class="w-full flex flex-col items-center mt-2">
  <img src="./images/image1.png"
    class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">
    초기 Brinkman (좌): 압력 강하 없음 ❌ /
    Darcy 표면항 (우): 압력 점프 재현 ✅
  </div>
</div>

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 11-1. 막 변위 결과

<div class="text-xl opacity-80 mb-6">압력 하중에 의한 탄성막 변형</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 막 변위 분포

<div class="w-full flex flex-col items-center mt-2">
  <img src="./images/image1.png"
    class="w-full max-h-[180px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">
    시간에 따른 막 변위 $w(r,t)$ — 중심부 최대 변위 발생
  </div>
</div>

</div>

<div v-click>

### 🔹 ALE 메시 변형

<div class="w-full flex flex-col items-center mt-2">
  <img src="./images/image1.png"
    class="w-full max-h-[180px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">
    막 변위에 따른 유체 메시 변형 (조화 확장)
  </div>
</div>

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 12. 수치적 한계 및 고찰

<div class="text-xl opacity-80 mb-6">알고 사용하는 근사들</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 한계 1: CG1 압력 공간
```python
Q_f = functionspace(msh, belem("Lagrange", ..., 1))
# CG1: 연속 공간 → 진짜 불연속 점프 표현 불가
```

- CG1은 메시 전체에서 압력이 연속적으로 이어져야 함
- 진짜 압력 불연속 표현: DG 압력 공간 필요
- 현재: 막 인근 요소에 걸쳐 **스미어된(Smeared) 압력 점프** 발생

### 🔹 한계 2: 1차 명시적 대류
```python
# Forward Euler 대류 (1차 정확도)
rho_v * inner(
    dot(u_n - w_msh, nabla_grad(u_n)), vf)
# 엄밀하게는 Adams-Bashforth 2차가 필요
# but: u_{n-1} 저장 필요 → 구현 복잡도 증가
```

고유속에서 수치적 불안정 가능성 존재

</div>

<div v-click>

### 🔹 한계 3: 명시적 결합 오차

스태거드 결합은 $O(\Delta t)$ 결합 오차 발생:

- NS → 막 변위 → ALE 순으로 차례로 풀기
- 엄밀하게는 Picard 내부 반복 필요
- 현재: $\Delta t$를 충분히 작게 하여 보완

### 🔹 한계 요약 및 대응

| 한계 | 현재 근사 | 영향 | 대응 |
|---|---|---|---|
| CG1 압력 | 스미어된 점프 | 정량 오차 | DG 공간 |
| 1차 대류 | Forward Euler | 안정성 | AB 2차 |
| 명시적 결합 | $O(\Delta t)$ 오차 | 결합 정확도 | Picard |
| 법선 샘플링 | 최소자승 근사 | 샘플링 오차 | dS 약형 |

<div class="mt-2 p-2 bg-blue-900 bg-opacity-20 rounded border border-blue-500 border-opacity-30 text-xs">
✅ 위 한계들은 알고 사용하는 근사이며 논문에 명시합니다.
</div>

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 13. 향후 연구 방향

<div class="text-xl opacity-80 mb-6">케이크층 모델링 및 개선 방향</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm">

<div v-click>

### 🔹 케이크층(Cake Layer) LVPP 적용

MBR 공정에서 슬러리 입자가 막 표면에 쌓여 케이크층을 형성합니다.

$$\frac{\partial c}{\partial t} = \underbrace{k_{dep} \cdot J_n}_{\text{퇴적}} - \underbrace{k_{ero} \cdot \tau_w}_{\text{침식}}, \quad c \geq 0$$

- $c$: 케이크층 두께
- $J_n = u \cdot n$: 막 법선 방향 유량 (퇴적 구동력)
- $\tau_w$: 전단 응력 (침식 구동력)
- $c \geq 0$: 자연스러운 부등식 제약 → **LVPP의 정확한 적용 대상**
```
막(Mesh로 표현) + 케이크층(LVPP로 표현)

→ 두 방법의 역할 분담이 물리적으로 올바름
→ 막: 위치가 명확 → 경계 조건
→ 케이크층: 위치 미지, 성장 → LVPP
```

</div>

<div v-click>

### 🔹 수치 정확도 개선

| 개선 항목 | 현재 | 목표 |
|---|---|---|
| 압력 공간 | CG1 | DG0/DG1 |
| 대류 이산화 | Forward Euler | Adams-Bashforth 2차 |
| 결합 방식 | 명시적 | Picard 반복 |

### 🔹 현재 시뮬레이션의 의미
```
BioWin (상용, 거시):
  공정 전체 성능 예측
         ↕ 상호 보완
FEniCSx (본 연구, 미시):
  막 내부 국소 압력·속도·변위
  케이크층 성장 동역학 (향후)
```

향후 공정 장치 설계 및 운전 최적화에 활용될 수 있는 **수치해석적 근거**를 제시합니다.

</div>
</div>
---
theme: seriph
class: text-center
highlighter: shiki
---

# 감사합니다

<div class="text-xl opacity-80 mb-8">Questions & Discussion</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm mt-4">

<div>

### 참고 문헌

- Donea, J. et al. (2004). Arbitrary Lagrangian-Eulerian Methods. *Encyclopedia of Computational Mechanics*.
- Keith, B. and Surowiec, T.M. (2024). Proximal Galerkin: A Structure-Preserving Finite Element Method for Pointwise Bound Constraints. *Found Comput Math*.

</div>

<div>

### 사용 소프트웨어

| 소프트웨어 | 용도 |
|---|---|
| FEniCSx v0.10.0 | 유한요소 계산 |
| Gmsh | 메시 생성 |
| ParaView | 결과 시각화 |
| PETSc / MPI | 선형 대수 / 병렬 |
| dolfinx | FEM 프레임워크 |

</div>

</div>

<div class="absolute bottom-6 right-6 text-xl">
  윤현준 / 박기성 / 박찬서 / 이민용
</div>
