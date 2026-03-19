# 나비에-스토크스 방정식 (IPCS 기법) 약형 변환 정리

비압축성 유체의 지배 방정식인 나비에-스토크스 방정식을 컴퓨터(FEniCSx)로 풀기 위해 IPCS 분할법으로 이산화하고 약형(Weak Form)으로 변환하는 과정입니다.

## 1. 지배 방정식 (Navier-Stokes)
유체의 운동량 보존과 비압축성을 나타냅니다.

* **운동량 보존 및 연속 방정식:**
$$\rho \left(\frac{\partial u}{\partial t} + u \cdot \nabla u\right) = \nabla \cdot \sigma(u, p) + f$$
$$\nabla \cdot u = 0$$
* **응력 텐서 (Stress Tensor):**
$$\sigma(u, p) = 2\mu\epsilon(u) - pI$$
$$\epsilon(u) = \frac{1}{2}(\nabla u + (\nabla u)^T)$$

---

## 2. Step 1: 임시 속도 $u^*$ 구하기
시간 이산화와 대류항 선형화를 거쳐 임시 속도를 구하는 단계입니다.

* **시간 이산화 (Crank-Nicolson):** 시간 미분 $\frac{\partial u}{\partial t}$를 $\frac{u^*-u^n}{\Delta t}$로 쪼개고, 점성항에 현재와 과거의 평균을 적용합니다.
* **대류항 선형화:** 비선형항 $u \cdot \nabla u$를 풀기 위해, 앞의 $u$는 과거 데이터를 이용한 Adams-Bashforth 외삽($u_{AB}$)으로, 뒤의 $\nabla u$는 Crank-Nicolson 미분($\nabla u_{CN}$)으로 처리합니다.
* **부분 적분 (약형 변환):** 양변에 테스트 함수 $v$를 곱하고 적분하여 2계 미분을 1계 미분으로 낮춥니다.

**최종 수식 ($F_1=0$):**
$$F_1 = \int \rho \frac{u^* - u^n}{\Delta t} \cdot v + \int \rho (u_{AB} \cdot \nabla u_{CN}) \cdot v + \int \mu \nabla u_{CN} : \nabla v - \int p^n (\nabla \cdot v) + \int f \cdot v = 0$$

---

## 3. 다공성 막 투과 모델링 (Darcy's Law 적용)

얇은 반투막(Membrane)을 통과하는 유동은 수만 개의 미세 기공을 3D 메쉬로 직접 묘사하는 것이 불가능합니다. 따라서 막을 2D 내부 경계면($\Gamma_m$)으로 정의하고, **Darcy의 법칙(Darcy's Law)**을 수학적인 표면 저항(Surface Resistance)으로 변환하여 적용합니다.

### 4. 물리적 배경 (Darcy 점프 조건)
다공성 매질을 통과하는 유속은 막 양단의 압력차($[[p]]$)에 비례하고 유체의 점성에 반비례합니다. 이를 막 표면에 수직으로 통과하는 법선 속도($u \cdot n$)에 대한 저항력으로 정리하면 다음과 같습니다.

$$[[p]] = \frac{\mu}{\kappa} (u \cdot n)$$
- $\kappa$: 막의 투과율 (Permeability)
- $\mu$: 유체의 점성 (Viscosity)
- $[[p]]$: 막을 기점으로 한 압력 강하 ($p^- - p^+$)

### 5. 약형 (Weak Form) 변환
이 압력 강하 조건은 유체가 막을 통과할 때 속도에 비례하여 받는 '마찰 저항력(Traction force)'으로 작용합니다. 따라서 Step 1의 운동량 방정식(임시 속도 구하기) 약형에서, 막 표면($\Gamma_m$)에 대한 면적분 항으로 추가됩니다.

양변에 테스트 함수 $v$를 내적하고 적분하면 다음과 같은 표면 저항항이 도출됩니다.
$$\int_{\Gamma_m} \frac{\mu}{\kappa} (u \cdot n) (v \cdot n) dS$$

### 6. FEniCSx 코드 구현
내부 경계면(Interior facet)에서는 유체가 통과할 때 면을 공유하는 양쪽 셀(Cell)의 물리량이 충돌할 수 있습니다. 수치적 안정성을 위해 속도의 평균값(`avg`)과 한쪽 방향의 법선 벡터(`n("+")`)를 취하여 계산합니다.

---

## 5. 움직이는 격자 기법 (ALE: Arbitrary Lagrangian-Eulerian)

유체-구조 연성(FSI) 문제에서 반투막과 같은 고체 구조물이 유압에 의해 변형되면, 유체가 흐르는 공간 자체의 형태가 변하게 됩니다. 이를 해결하기 위해 **ALE(임의 라그랑지안-오일러리안)** 기법을 도입합니다.

---

### 1. 물리적 배경: 왜 ALE를 쓰는가?
유체역학과 고체역학은 격자(Mesh)를 바라보는 관점이 다릅니다.
* **Eulerian (오일러 관점, 유체):** 격자는 공간에 고정되어 있고, 유체가 그 격자를 통과하여 지나갑니다. (다리 위에서 강물을 바라보는 관점)
* **Lagrangian (라그랑주 관점, 고체):** 격자가 물질(입자)과 함께 움직입니다. 고체의 변형을 추적하기 좋지만, 유체에 적용하면 격자가 심하게 꼬여서(Tangling) 계산이 파탄 납니다.
* **ALE (절충형):** 구조물의 경계면에서는 격자가 고체와 함께 움직이고(Lagrangian), 유체 내부에서는 격자가 꼬이지 않도록 부드럽게 재배치되며(Eulerian), 그 사이의 공간은 '임의의 속도'로 움직입니다.

### 2. 수식 변환: 움직이는 메쉬에서의 나비에-스토크스
격자 자체가 $w_{mesh}$라는 속도로 움직이고 있기 때문에, 유체가 느끼는 '상대적인 대류 속도'가 달라집니다. 따라서 기존 나비에-스토크스 방정식의 대류항(Convection term)이 다음과 같이 수정됩니다.

**기존 대류항:**
$$\rho (u \cdot \nabla u)$$
**ALE 대류항 (메쉬 속도 $w_{mesh}$ 반영):**
$$\rho ((u - w_{mesh}) \cdot \nabla u)$$

- $u$: 유체의 실제 속도
- $w_{mesh}$: 격자(Mesh)가 이동하는 속도
- $(u - w_{mesh})$: 움직이는 격자 위에서 유체가 움직이는 상대 속도

### 3. 메쉬 변위의 조화 확장 (Harmonic Extension)
경계면(막)이 움직일 때, 내부 유체 격자들이 서로 겹치거나 꼬이지 않도록 변위($d$)를 공간 전체로 부드럽게 퍼뜨려야 합니다. 이를 위해 공간 전체에 대해 라플라스 방정식(Laplace's equation)을 풉니다.

$$-\nabla \cdot (k \nabla d) = 0$$
- $d$: 메쉬의 변위 (Mesh displacement)
- $k$: 메쉬 강성 (Mesh stiffness). 막에 가까울수록 크게 설정하여 좁은 격자의 찌그러짐을 방지합니다.

### 4. FEniCSx 코드 구현
ALE 기법은 유체 운동량 방정식의 대류항 수정과, 메쉬 변위를 구하는 Poisson 방정식 두 가지로 구현됩니다.

```python
# 1. ALE 대류항 적용 (u_n - w_msh)
F1 = rho_v / dt * inner(u_t - u_n, vf) * dxf \
   + rho_v * inner(dot(u_n - w_msh, nabla_grad(u_n)), vf) * dxf \
   + inner(sig_(U, p_n), eps_(vf)) * dxf

# 2. 메쉬 꼬임 방지를 위한 가변 강성(Variable Stiffness) 조화 확장
dist_to_mem = abs(x_f[0] - x_m)
mesh_stiffness = 1.0 / (dist_to_mem + 0.05)**2 

a_ale_ufl = form(mesh_stiffness * inner(nabla_grad(d_t3), nabla_grad(e3)) * dxf)
L_ale_ufl = form(inner(Constant(msh, PETSc.ScalarType((0., 0., 0.))), e3) * dxf)