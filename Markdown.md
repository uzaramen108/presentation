# 유체-막 상호작용 시뮬레이션
## Navier-Stokes → Brinkman → Darcy + ALE 로의 발전

---

## 목차
1. [시뮬레이션 개요](#1-시뮬레이션-개요)
2. [나비에-스토크스 약형 유도](#2-나비에-스토크스-약형-유도-ipcs-기법)
3. [막 투과 모델: Brinkman 저항항](#3-막-투과-모델-brinkman-저항항)
4. [Brinkman 방식의 수학적 한계](#4-brinkman-방식의-수학적-한계)
5. [올바른 해법: Navier-Stokes + Darcy](#5-올바른-해법-navier-stokes--darcy)
6. [ALE 기법](#6-ale-기법-arbitrary-lagrangian-eulerian)

---

## 1. 시뮬레이션 개요

### 물리적 설정

원형 단면의 3D 파이프 내부를 유체가 흐르다가 중간에 위치한 **반투막(semi-permeable membrane)**과 만나는 문제입니다.

| 파라미터 | 값 |
|---|---|
| 파이프 길이 $L$ | 5 m |
| 파이프 반지름 $R$ | 0.5 m |
| 막 위치 $x_m$ | 3.75 m |
| 유체 밀도 $\rho$ | 1.0 kg/m³ |
| 유체 점성 $\mu$ | 0.01 Pa·s |
| 막 투과율 $\kappa$ | 5×10⁻³ m² |
| 막 장력 $T$ | 25 N/m |

입구 속도는 Poiseuille 포물선 프로파일에 사인파 시간 변화를 적용합니다:

$$u_{inlet}(r,t) = u_{max} \sin\!\left(\frac{\pi t}{2T}\right)\left(1 - \frac{r^2}{R^2}\right)$$

### 계산 목표

나비에-스토크스 방정식으로 매 타임스텝마다 다음을 계산합니다.

- **속도장** $u(x,t)$: 막 전후 유체의 흐름 패턴
- **압력장** $p(x,t)$: 막 전후 압력 분포
- **막 변위** $w(x,t)$: 압력 하중에 의한 막의 휨

### 수치 기법: IPCS 분할법

나비에-스토크스를 직접 풀면 속도-압력 연성(Coupling)으로 인해 계산 비용이 극도로 높아집니다. 이를 **3단계로 분리**하여 푸는 IPCS(Incremental Pressure Correction Scheme)를 사용합니다.

$$\text{Step 1: } u^* \;\rightarrow\; \text{Step 2: } p^{n+1} \;\rightarrow\; \text{Step 3: } u^{n+1}$$

---

## 2. 나비에-스토크스 약형 유도 (IPCS 기법)

### 2.1 지배 방정식 (Strong Form)

비압축성 나비에-스토크스 운동량 방정식과 연속 방정식:

$$\rho \frac{\partial u}{\partial t} + \rho(u \cdot \nabla u) = \nabla \cdot \sigma(u,p) + f$$

$$\nabla \cdot u = 0$$

응력 텐서와 변형률 텐서:

$$\sigma(u,p) = 2\mu\epsilon(u) - pI, \qquad \epsilon(u) = \frac{1}{2}\left(\nabla u + (\nabla u)^T\right)$$

---

### 2.2 Step 1: 임시 속도 $u^*$ 구하기

#### ① 시간 이산화 (Crank-Nicolson)

시간 미분 $\frac{\partial u}{\partial t}$를 $\frac{u^*-u^n}{\Delta t}$로 치환하고, 점성항에 현재와 과거의 평균을 적용합니다.

#### ② 대류항 선형화

비선형항 $u \cdot \nabla u$를 풀기 위해:

- **앞부분 — Adams-Bashforth 외삽**: $u_{AB} = \dfrac{3}{2}u^n - \dfrac{1}{2}u^{n-1}$
- **뒷부분 — Crank-Nicolson 평균**: $u_{CN} = \dfrac{u^*+u^n}{2}$

이산화된 운동량 방정식:

$$\rho\frac{u^* - u^n}{\Delta t} + \rho(u_{AB}\cdot\nabla u_{CN}) - \nabla\cdot\!\left(2\mu\epsilon(u_{CN})\right) + \nabla p^n = f$$

#### ③ 약형 변환 (Weak Form)

테스트 함수 $v$를 곱하고 $\Omega$에서 적분한 뒤, **부분 적분(Integration by Parts)**으로 미분 차수를 낮춥니다.

$$\int_\Omega \nabla\cdot(2\mu\epsilon)\cdot v \;\xrightarrow{\text{부분 적분}}\; \int_\Omega 2\mu\,\epsilon(u_{CN}):\epsilon(v)$$

$$\int_\Omega \nabla p \cdot v \;\xrightarrow{\text{부분 적분}}\; -\int_\Omega p(\nabla\cdot v)$$

#### ④ 최종 약형 $F_1 = 0$

$$\underbrace{\int_\Omega \rho\frac{u^*-u^n}{\Delta t}\cdot v}_{\text{시간 미분}} + \underbrace{\int_\Omega \rho(u_{AB}\cdot\nabla u_{CN})\cdot v}_{\text{대류항}} + \underbrace{\int_\Omega 2\mu\,\epsilon(u_{CN}):\epsilon(v)}_{\text{점성항}} - \underbrace{\int_\Omega p^n(\nabla\cdot v)}_{\text{압력항}} + \underbrace{\int_\Omega f\cdot v}_{\text{외력}} = 0$$

---

### 2.3 Step 2: 압력 보정 (Pressure Correction)

진짜 속도 $u^{n+1}$과 임시 속도 $u^*$의 차이는 압력 증분 $\phi = p^{n+1} - p^n$에서 발생합니다:

$$\rho\frac{u^{n+1}-u^*}{\Delta t} = -\nabla\phi$$

양변에 $\nabla\cdot$을 취하고 $\nabla\cdot u^{n+1}=0$을 대입하면 **압력 포아송 방정식**이 됩니다:

$$\nabla^2\phi = \frac{\rho}{\Delta t}\nabla\cdot u^*$$

테스트 함수 $q$를 곱하고 부분 적분하면 최종 약형:

$$\int_\Omega \nabla\phi\cdot\nabla q = -\int_\Omega \frac{\rho}{\Delta t}(\nabla\cdot u^*)\,q$$

---

### 2.4 Step 3: 속도 보정 (Velocity Projection)

Step 2에서 구한 $\phi$로 최종 속도를 보정합니다:

$$u^{n+1} = u^* - \frac{\Delta t}{\rho}\nabla\phi$$

약형으로 변환하면:

$$\int_\Omega \rho\,u^{n+1}\cdot v = \int_\Omega \rho\,u^*\cdot v - \int_\Omega \Delta t\,\nabla\phi\cdot v$$

이 단계를 거친 $u^{n+1}$은 비압축성 조건 $\nabla\cdot u = 0$을 정확히 만족합니다. 압력도 다음과 같이 갱신됩니다:

$$p^{n+1} = p^n + \phi$$

---

## 3. 막 투과 모델: Brinkman 저항항

### 3.1 Brinkman 방정식

다공성 매질 내부에서의 유동을 기술하는 방정식입니다. 나비에-스토크스에 **부피 저항항** $\dfrac{\mu}{\kappa}u$를 추가한 형태입니다:

$$\rho\frac{\partial u}{\partial t} + \rho(u\cdot\nabla u) = -\nabla p + \mu\nabla^2 u - \underbrace{\frac{\mu}{\kappa}u}_{\text{Brinkman 저항}}$$

- $\kappa$: 투과율 (Permeability)
- $\dfrac{\mu}{\kappa}$: 물리적 저항 계수 (투과율이 작을수록 저항 증가)

### 3.2 코드에서의 적용 방식

막 위치($|x - x_m| < \delta$)에서만 저항이 활성화되는 조건부 함수로 구현합니다:

$$F_1 \;\mathrel{+}=\; \int_\Omega \frac{\mu}{\kappa}\cdot\mathbf{1}_{\Gamma_m^\delta}\cdot u\cdot v\;dx$$

$\mathbf{1}_{\Gamma_m^\delta}$는 막 근방 영역에서만 1, 나머지에서 0인 지시함수입니다.

### 3.3 저항 계수 설정

| $\kappa \to \infty$ | $\kappa \to 0$ |
|---|---|
| 저항 $\to 0$ (자유 통과) | 저항 $\to \infty$ (완전 차단) |

초기 구현에서는 $R_m = 50$을 임의로 설정했으나, 이후 물리량 기반으로 $R_m = \mu/\kappa$로 수정하였습니다.
<video src="./video/Brinkman-velocity.mp4" controls></video>

---

## 4. Brinkman 방식의 수학적 한계

### 4.1 물리적으로 기대되는 현상

반투막을 유체가 통과할 때, Darcy의 법칙에 의해 막 전후로 **압력 강하(Pressure Jump)**가 반드시 발생해야 합니다:

$$\Delta p = p^- - p^+ = \frac{\mu}{\kappa}(u\cdot n)\cdot d_{mem}$$

- 막 전(상류): 압력 $p^-$ (높음)
- 막 후(하류): 압력 $p^+$ (낮음)

### 4.2 ParaView 시각화 결과

Brinkman 방식으로 시뮬레이션한 결과를 ParaView로 확인하면, **막 전후 압력 분포가 거의 동일**하게 나타납니다. 압력이 막을 기점으로 불연속적으로 떨어지는 현상(압력 점프 $[\![p]\!]$)이 관찰되지 않습니다.
<video src="./video/Brinkman-pressure.mp4" controls></video>

### 4.3 수학적 원인

Brinkman 항은 **부피 적분(Volume integral)**으로 추가됩니다:

$$\int_{\Omega} \frac{\mu}{\kappa}\cdot\mathbf{1}_{\Gamma_m^\delta}\cdot u\cdot v\;dx \qquad \leftarrow \text{3D 영역 전체에 분산}$$

그러나 실제 막에서의 압력 점프는 **표면 적분(Surface integral)**이어야 합니다:

$$\int_{\Gamma_m} [\![p]\!]\cdot v\cdot n\;dS \qquad \leftarrow \text{2D 경계면에 집중}$$

> **결론**: Brinkman은 막 근처에서 유체를 느리게 만드는 효과는 있지만, 막 양면의 압력 불연속(Pressure Jump)을 물리적으로 올바르게 표현하지 못합니다.

---

## 5. 올바른 해법: Navier-Stokes + Darcy

### 5.1 Darcy 점프 조건 (Surface Resistance)

막을 3D 부피가 아닌 **2D 내부 경계면** $\Gamma_m$으로 정의하면, 막 통과 유동을 Darcy 법칙으로 정확히 기술할 수 있습니다:

$$[\![p]\!] = \frac{\mu}{\kappa}(u\cdot n), \qquad [\![p]\!] = p^- - p^+$$

이를 NS 운동량 방정식 약형에 **표면 적분항**으로 추가합니다:

$$F_1 \;\mathrel{+}=\; \int_{\Gamma_m} \frac{\mu}{\kappa}\;\text{avg}(u)\cdot n^+\;\text{avg}(v)\cdot n^+\;dS$$

- `avg`: 막 양쪽 셀의 평균값 (수치 안정성 확보)
- $n^+$: 막 한쪽 방향의 법선 벡터
- 계수 $\mu/\kappa$: 투과율에서 유도된 물리량 (임의 상수가 아님)

### 5.2 Brinkman vs Darcy 비교

| | Brinkman | Darcy (표면항) |
|---|:---:|:---:|
| 적분 형태 | 부피 $dx$ | 표면 $dS$ |
| 압력 점프 표현 | ❌ | ✅ |
| 물리적 근거 | 근사 | 엄밀 |
| 구현 난이도 | 낮음 | 높음 |
| 메시 요구사항 | 단일 메시 | 내부 경계면 필요 |

### 5.3 한계: 막 변형을 다루려면 ALE가 필요

막이 압력을 받아 **변형(변위)**되면, 유체 도메인 $\Omega_f(t)$ 자체가 시간에 따라 변합니다. 고정된 메시 위에서 NS를 풀면 변형된 형상을 인식하지 못하므로, **움직이는 경계를 다루는 ALE 기법**이 추가로 필요합니다.

---

## 6. ALE 기법 (Arbitrary Lagrangian-Eulerian)

### 6.1 왜 ALE가 필요한가?

유체역학과 고체역학은 격자(Mesh)를 바라보는 관점이 서로 다릅니다.

| 관점 | 설명 | 한계 |
|---|---|---|
| **Eulerian** (유체역학) | 격자 고정, 유체가 격자를 통과 | 움직이는 경계 추적 불가 |
| **Lagrangian** (고체역학) | 격자가 물질과 함께 이동 | 유체에 적용 시 격자 뒤틀림(Tangling) |
| **ALE** (절충형) | 경계는 Lagrangian, 내부는 Eulerian | — |

ALE 기법에서:
- **구조물 경계면(막)**에서는 격자가 고체와 함께 이동합니다 (Lagrangian).
- **유체 내부**에서는 격자가 꼬이지 않도록 부드럽게 재배치됩니다 (Eulerian).

### 6.2 ALE 나비에-스토크스

격자 자체가 $w_{mesh}$의 속도로 움직이므로, 유체가 느끼는 **상대적 대류 속도**가 달라집니다.

**기존 대류항:**
$$\rho(u\cdot\nabla u)$$

**ALE 대류항 (메시 속도 반영):**
$$\rho\underbrace{(u - w_{mesh})}_{\text{상대 속도}}\cdot\nabla u$$

전체 ALE-NS 약형 ($F_1 = 0$):

$$\int_\Omega \rho\frac{u^*-u^n}{\Delta t}\cdot v + \int_\Omega \rho(u^n - w_{mesh})\cdot\nabla u^n\cdot v + \int_\Omega \sigma:\epsilon(v)\;dx + \int_{\Gamma_m}\frac{\mu}{\kappa}\,\text{avg}(u)\cdot n\;\text{avg}(v)\cdot n\;dS = 0$$

### 6.3 조화 확장 (Harmonic Extension)

막이 변위 $w$만큼 움직이면, 내부 격자들이 겹치지 않도록 변위를 **라플라스 방정식**으로 도메인 전체에 부드럽게 전파합니다:

$$-\nabla\cdot(k\nabla d) = 0 \quad \text{in } \Omega_f$$

경계 조건:

$$d = \Delta w\,\hat{e}_x \;\text{ on } \Gamma_m \qquad d = 0 \;\text{ on } \partial\Omega_f$$

- $d$: 메시 변위 (Mesh displacement)
- $\Delta w$: 막의 증분 변위 (이번 스텝에서 새로 발생한 변위)
- $k$: **가변 강성(Variable Stiffness)** — 막에 가까울수록 크게 설정하여 좁은 격자의 뒤틀림을 방지합니다:

$$k = \frac{1}{(|x - x_m| + \delta)^2}$$

### 6.4 전체 타임스텝 결합 구조

매 타임스텝마다 다음 순서로 연성 계산이 수행됩니다:
```
① ALE-NS 풀기 (u*, p → u^{n+1}, p^{n+1})
         ↓
② 막 양단 압력 추출: [[p]] = p⁻ - p⁺
         ↓
③ 막 변위 (Poisson): -T Δ_Γ w = [[p]]
         ↓
④ 메시 변위 (조화 확장): -∇·(k∇d) = 0
         ↓
⑤ 메시 좌표 갱신: x^{n+1} = x^n + d
         ↓
⑥ bb_tree 갱신 후 다음 스텝으로
```

이 구조를 통해 유체 속도·압력, 막 변위, 메시 형상이 **물리적으로 일관되게 결합**되어 시뮬레이션됩니다.

---

## 참고 문헌

- Beavers, G.S. & Joseph, D.D. (1967). Boundary conditions at a naturally permeable wall. *Journal of Fluid Mechanics*, 30(1), 197–207.
- Donea, J. et al. (2004). Arbitrary Lagrangian-Eulerian Methods. *Encyclopedia of Computational Mechanics*.
- Langtangen, H.P. & Logg, A. (2016). *Solving PDEs in Python: The FEniCS Tutorial*. Springer.
- Chorin, A.J. (1968). Numerical solution of the Navier-Stokes equations. *Mathematics of Computation*, 22(104), 745–762.