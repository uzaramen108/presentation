---
theme: seriph
background: ./images/white_background.png
title: LVPP 알고리즘을 활용한 장애물 문제 분석
info: |
  ## Proximal Galerkin & Obstacle Problem
  3월 17일 발표 자료
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
duration: 20min
---

# LVPP 알고리즘을 활용한
# 장애물 문제(Obstacle Problem) 분석

<div class="text-xl opacity-80 mb-8">
  Proximal Galerkin 기반의 구조 보존형 유한요소법
</div>

<div class="abs-br m-6 text-xl">
  윤현준 / 박기성 / 박찬서 / 이민용
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 1. LVPP의 도입 배경

<div class="text-xl opacity-80 mb-8">
  부등식 제약 조건이 유발하는 Nonsmoothness
</div>

<div grid="~ cols-2 gap-8" class="text-left">

<div v-click>

### 1. 제약이 포함된 변분 문제
과학과 공학의 다양한 분야에서는 부등식 제약 조건(Inequality constraints)을 포함하는 변분 문제가 빈번하게 발생합니다.

<div v-click class="mt-4 text-sm opacity-80">

- 농도가 음수가 될 수 없음 $(C\ge 0)$
- 몰 분율은 0부터 1 사이의 값을 지님 $(0\leq x_{i}\leq 1, \sum x_{i} = 1)$
- 특정 stress 이상에서 소성 변형이 일어남 $(\nabla u \leq \phi)$

</div>

</div>

<div v-click>

### 2. 기존 PDE 솔버의 한계
부등식 제약은 본질적으로 **Nonsmoothness**를 유발하여 일반적인 편미분 방정식(PDE) 솔버로 풀기 매우 까다롭습니다.

<div class="mt-4 text-sm opacity-80">

- **한계** : 최적성 조건이 연속 차원에서 충분히 매끄럽지 않음
- **해결** : 이를 효율적으로 풀기 위해 **LVPP (Latent Variable Proximal Point)** 프레임워크 도입

</div>

<div class="w-full flex justify-center mt-6">
  <img src="./images/image_1.png" alt="Obstacle Problem Membrane" class="w-full max-h-[160px] object-contain rounded shadow-sm" />
</div>

</div>
</div>---
theme: seriph
class: text-center
highlighter: shiki

---
theme: seriph
class: text-center
highlighter: shiki
---

# 2. 기존 수치해석의 한계

<div class="text-xl opacity-80 mb-8">
  메시 의존성과 비벡터공간 문제
</div>

<div grid="~ cols-2 gap-8" class="text-left">

<div v-click>

### 1. 메시 의존성 (Mesh-dependence)
부등식 제약 문제는 연속 차원에서 충분히 매끄럽지 않습니다. 

<div v-click class="mt-4 text-sm opacity-80">

- **한계점** : 연속 수준에서 Semismooth Newton 방법 등을 바로 적용하기 어려움.
- **성능 저하** : 이산화 과정에서 격자(Mesh)를 세밀하게 나누거나 고차(High-order) 다항식을 사용할수록, 선형 시스템을 푸는 반복 연산 횟수가 기하급수적으로 폭증함.

</div>

</div>

<div v-click>

### 2. 비벡터공간 (Non-vector space)
부등식 제약 때문에 해를 만족하는 집합(Feasible set $K$)은 벡터 공간을 이루지 못합니다.

<div class="mt-4 text-sm opacity-80">

- **제약 이탈** : 해를 선형 결합하거나 상수를 곱하면 제약 조건<br>($u \ge \phi$)을 벗어날 수 있음.
- **적용 불가** : 따라서 표준 유한요소법(FEM)의 기저 함수 선형 결합이나, Newton 방법의 가산적 업데이트(Additive update) 방식을 직접 적용할 수 없음.

</div>

<div class="w-full flex justify-center mt-6">
  <img src="./images/image_2.png" alt="Mesh Dependence Chart" class="w-full max-h-[160px] object-contain rounded shadow-sm" />
</div>

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 3. LVPP의 핵심 아이디어

<div class="text-xl opacity-80 mb-8">
  잠재 변수와 매끄러운 좌표 변환
</div>

<div grid="~ cols-2 gap-8" class="text-left">

<div>

### 1. 매끄러운 좌표 변환
제약 조건 집합 $K$를 직접 다루는 대신, 제한이 없는 무한 차원의 잠재 함수 공간(Latent function space) $\psi$를 도입합니다.

<div class="mt-4 text-sm opacity-80">

- **Legendre function** : 르장드르 함수 $R(u)$를 이용하여 제약 공간과 잠재 공간을 연결하는 매끄러운 좌표 변환을 생성합니다.
- **변수 정의** : $\psi = \nabla R(u)$
- **예시($u\ge \phi$의 경우)** : $R(u) = (u-\phi)\ln(u-\phi)-(u-\phi)$
</div>

</div>

<div v-click>

### 2. 구조 보존 (Structure-preserving)
잠재 변수에서 물리적 변위 $u$로 돌아오는 역변환 식을 적용합니다.

$$
u = \nabla R^*(\psi) = \phi + e^\psi
$$

<div class="mt-4 text-sm opacity-80">

- **자동 제약 만족** : 지수 함수($e^\psi > 0$)가 포함되어 있어, $\psi$가 어떤 실수 값을 갖더라도 $u$는 무조건 장애물 $\phi$보다 큽니다.
- **알고리즘적 이점** : 제약이 물리적 벽처럼 보존되므로, 표준 FEM 기저 함수의 선형 결합 및 가산적 업데이트(Additive update)가 가능해집니다.

</div>

<div class="w-full flex justify-center mt-6">
  <img src="./images/image_3.png" alt="Coordinate Transformation" class="w-full max-h-[140px] object-contain rounded shadow-sm" />
</div>


</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 4. LVPP의 주요 특징
<div grid="~ cols-2 gap-8" class="text-left">

<div class="mt-4 text-sm opacity-80">

### 👍 장점 (Pros)
* **무한 차원 정식화**: 격자 크기(Mesh size)나 이산화 차수(Polynomial order)에 관계없이 성능이 유지되는 메시 및 차수 독립성을 가집니다.
* **투영(Projection) 불필요**: 이산 수준에서 점별 제약(Pointwise constraints)을 강제하는 단순한 메커니즘을 제공하여, 엄격하게 제약조건을 준수합니다.
* **구현 용이성 및 호환성**: 특수한 이산화 기법 없이, 고차 수치 방법 및 FEniCSx와 같은 표준 유한요소 라이브러리와 완벽하게 호환됩니다.
* **강건한 수치 성능**: 파라미터를 무한대로 보내는 과정에 의존하지 않고 안정적으로 수렴합니다.

</div>

<div v-click>
<div class="mt-4 text-sm opacity-80">

### 👎 단점 및 고려사항 (Trade-offs)
* **근위 매개변수 $\alpha_k$ 설정 딜레마**:
  * 최적화 가속을 위해 $\alpha_k$를 지수적으로 너무 공격적으로 키우면, 내부 비선형 PDE 서브프라블럼이 악조건(Ill-conditioned) 상태가 되어 Newton 솔버 수렴이 어려워집니다.
  * 반대로 너무 천천히 키우면 전체 LVPP 반복 횟수가 늘어나는 트레이드오프가 존재합니다.
* **연산 비용 증가**: 선형 시스템을 한 번 푸는 것으로 끝나지 않고, 매 스텝마다 비선형 대수 방정식 시스템을 연속적으로 반복해서 풀어야 하므로 초기 세팅 및 단계별 연산 비용이 높습니다.

</div>

<div class="w-full flex justify-center mt-6">
  <img src="./images/image_3.png" alt="Coordinate Transformation" class="w-full max-h-[140px] object-contain rounded shadow-sm" />
</div>


</div>
</div>


---
theme: seriph
class: text-center
highlighter: shiki
---

# 5. 장애물 문제(Obstacle Problem) 소개
<div class="text-xl opacity-80 mb-6">Dirichlet Energy 최소화 수식 모델링</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click >

### 🔹 원래의 문제 정의
장애물 문제는 다음의 디리클레 에너지 최소화 식으로 나타납니다 .
첫 번째 항은 막의 탄성 에너지, 두 번째 항은 외력을 의미합니다 .

$$
J(u) = \frac{1}{2} \int_\Omega \nabla u \cdot \nabla u \, dx - \int_\Omega f u \, dx
$$

막이 장애물 $\phi$를 통과할 수 없으므로 $u \ge \phi$라는 위치 제약조건이 발생하며, 최종 목표는 다음 제약 조건을 만족하는 해를 찾는 것입니다 .

$$
\min_{u \in K} J(u), \quad K = \{v \in H_0^1(\Omega) \mid v \ge \phi\}
$$

</div>
<div v-click >

### 🔹 조건 부연 설명
* **Feasible Set $K$**: 부등식 제약 조건을 만족하는 해의 집합을 의미합니다 .
* **경계 조건**: $\partial\Omega$ (경계면)에서 $u = 0$ 이고, 장애물 위에 놓인 상태를 가정합니다 .
* 제약 조건이 없다면 1차 변분을 통해 푸아송 방정식(Poisson equation)의 약형(Weak form)과 수학적으로 완전히 동일해집니다 .

</div>
</div>
</div>

<div v-click class="w-full flex justify-center mt-2">
  <img src="./images/image_4.png" alt="Dirichlet Energy Minimization" class="w-full max-h-[180px] object-contain rounded shadow-sm" />
</div>
---
theme: seriph
class: text-center
highlighter: shiki
---

# 5-1. 장애물 및 막의 수학적 정의
<div class="text-xl opacity-80 mb-6">벤치마크 모델의 도메인 및 장애물 수식</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click>

### 🔹 도메인 및 외력 조건
막이 존재하는 2차원 공간 $\Omega$와 외력 $f$를 다음과 같이 정의합니다 .

$$
\Omega = \{(x, y) \in \mathbb{R}^2 : 0 < r < 1\}, \quad r^2 = x^2 + y^2
$$
$$
f \equiv 0
$$

### 🔹 장애물 $\phi(x, y)$의 수식
막이 뚫을 수 없는 장애물은 중심부 반구 형태와 매끄럽게 이어지는 바깥 부분으로 정의됩니다 .

$$
\phi(x,y) = \begin{cases} \sqrt{1/4 - r^2} & r \le b, \\ d + b^2/d - br/d & r > b, \end{cases}
$$
※ $b = 9/20, \quad d = \sqrt{1/4 - b^2}$ 

</div>
<div v-click class="flex flex-col items-center">

### 🔹 장애물 단면 및 막의 거동
* 장애물의 중심부($r \le b$)는 위로 둥글게 솟아오른 구면을 이룹니다 .
* 탄성 막($u$)은 이 장애물 표면($\phi$)을 통과할 수 없으며 이는  $u \ge \phi$ 제약으로 나타납니다.

<div class="w-full flex flex-col items-center justify-center mt-6">
  <img src="./images/image_1.png" alt="Obstacle Cross Section" class="w-full max-h-[180px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-2">장애물 $\phi(x,y)$의 2D 단면 및 막의 접촉 형태</div>
</div>

</div>
</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 6-1 ~ 6-2. 르장드르 함수 및 잠재 변수 치환
<div class="text-xl opacity-80 mb-6">부등식 제약을 매끄러운 함수로 변환</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click >

### 🔹 섀넌 엔트로피 기반 변환
부등식 제약 조건 $u - \phi > 0$을 풀기 위해 섀넌 엔트로피(Shannon entropy) 기반의 르장드르 함수 $R(u)$를 정의합니다.

$$
R(u) = (u - \phi) \ln(u - \phi) - (u - \phi)
$$

### 🔹 잠재 변수 유도 및 역치환
위 식을 미분하여 제약이 없는 새로운 잠재 변수(latent constant)<br> $\psi$를 정의합니다. 이를 역으로 치환하면 물리적 해 $u$와 잠재 변수 $\psi$ 사이의 관계식이 도출됩니다.

$$
\psi = \nabla R(u) = \ln(u - \phi) \Rightarrow  u = \nabla R^*(\psi) = \phi + e^\psi
$$

**결과**: $\psi$가 모든 실수 값을 가질 수 있지만, 역치환된 변위 $u$는 항상 장애물보다 큰 값을 유지하게 됩니다.

</div>
<div v-click >

<div class="w-full flex flex-col items-center">
    <div class="text-xs opacity-70 mb-1">변환 전: 하한선이 존재하는 부등식 제약 공간 ($u \ge \phi$)</div>
    <img src="./images/image_1.png" alt="Before Legendre" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  </div>

  <div class="w-full flex flex-col items-center">
    <div class="text-xs opacity-70 mb-1">변환 후: 제약이 없는 무한한 실수 공간 ($\psi \in \mathbb{R}$)</div>
    <img src="./images/image_2.png" alt="After Legendre" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  </div>

</div>
</div>
</div>

<div v-click class="w-full flex justify-center mt-2">
  <img src="./images/image_4.png" alt="Dirichlet Energy Minimization" class="w-full max-h-[180px] object-contain rounded shadow-sm" />
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 6-3 ~ 6-4. Bregman 근위점 및 안장점 변환
<div class="text-xl opacity-80 mb-6">최적화 안정성 확보 및 라그랑지안 유도</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click >

### 🔹 Bregman Proximal Point Algorithm
원래의 최소화 문제를 **브레그만 발산(Bregman Divergence)** $D_R$을 이용해 변환합니다 . 이는 유클리드 거리 기반 방법을 보다 일반적인 기하 구조로 확장한 것입니다 .

$$
D_R(a,b) = R(a) - R(b) - \nabla R(b) \cdot (a-b)
$$

$$
u^k \in \arg \min_{u \in K} \left[ J(u) + \alpha_k^{-1} \int_\Omega D_R(u, u^{k-1}) \, dx \right]
$$

* **역할**: $D_R$은 거리와 유사한 역할을 하여, 이전 해 $u^{k-1}$에서 너무 멀리 가지 않도록 제한(Penalty)합니다 .
* **효과**: 각 반복(Iteration)마다 매개변수 $\alpha_k$를 조절함으로써 알고리즘의 **높은 안정성과 빠른 수렴**을 동시에 확보합니다 .

</div>
<div v-click >

### 🔹 안장점(Saddle Point) 문제로의 변환
LVPP 알고리즘의 핵심은 부등식 문제를 풀기 쉬운 **안장점 문제**로 바꾸는 것입니다 . 

* **Convex Conjugate**: 볼록 켤레 성질 <br>$R(u) = \sup_\psi (\langle \psi, u \rangle - R^*(\psi))$을 적용합니다 .
* **상수항 제거 및 치환**: 최적화에 영향을 주지 않는 상수항 $R(u^{k-1})$을 제거하고, $\psi$에 대해 변수 치환을 수행합니다 .

위 과정을 거치면 다음과 같은 라그랑지안 $L(u, \psi)$이 도출됩니다 .

$$
L(u, \psi) = J(u) + \alpha_k^{-1} \int_\Omega (\psi \cdot u - R^*(\psi + \nabla R(u^{k-1}))) \, dx
$$

최종적으로 $\min_u \max_\psi L(u, \psi)$ 구조에서 편미분($\frac{\partial L}{\partial \psi}=0, \frac{\partial L}{\partial u}=0$)을 통해 안장점 조건을 찾습니다 .

</div>
</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 6-5. 최종 약형식(Weak Form) 도출
<div class="text-xl opacity-80 mb-6">FEniCSx 구현을 위한 최종 연립방정식</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click >

### 🔹 에너지 균형 및 제약 조건 매핑
안장점의 편미분($\frac{\partial L}{\partial \psi} = 0$, $\frac{\partial L}{\partial u} = 0$) 결과를 정리하면 다음 두 식을 얻습니다 .

$$
u^k - \nabla R^*(\psi^k) = 0 \quad \text{--- (제약 조건 매핑)}
$$
$$
\alpha_k J'(u^k) + \psi^k = \psi^{k-1} \quad \text{--- (에너지 균형)}
$$

</div>
<div v-click >

### 🔹 최종 연립 약형식 (Weak Form)
위 식에 시험 함수 $v, w$를 곱하고 푸아송 약형식($J'(u^k)$)을 대입하여 연립 비선형 편미분 방정식의 약형식을 완성합니다 .

$$
(1) \quad \alpha_k(\nabla u^k, \nabla v) + (\psi^k, v) = \alpha_k(f, v) + (\psi^{k-1}, v)
$$
$$
(2) \quad (u^k, w) - (\phi + e^{\psi^k}, w) = 0
$$

결과적으로 다루기 힘들었던 부등식 최적화 문제가, FEniCSx로 구현하기 쉬운 **부드러운 안장점 문제**로 깔끔하게 변환되었습니다 .

</div>
</div>
</div>

---
layout: two-cols
layoutClass: gap-8
---

# 7-1. github 코드 실행 

::left::

<div class="text-sm mt-4">

### 🔹 예제 파일 구성
예제 파일은 다음과 같이 구성되어있으며 터미널(bash)에서 명령어로  실행합니다(docker).

### 🔹 각 파일 개요
* **generate_mesh_gmsh.py**: 장애물 및 막 메쉬 형성
* **compare_all.py**: 각 obstacle_-.py를 호출하여 계산된 결과를 비교.
* **obstacle_-.py**: 각 코드에 따른 솔버를 통하여 장애물 문제를 계산 및 저장.
</div>

::right::
<div v-click>
<div class="overflow-y-auto max-h-[450px] shadow-lg rounded-md border border-gray-200/20 text-sm mt-8">

```bash {all} twoslash
python3 generate_mesh_gmsh.py
python3 compare_all.py -P ./meshes/disk_1.xdmf -O coarse
python3 compare_all.py -P ./meshes/disk_2.xdmf -O medium
python3 compare_all.py -P ./meshes/disk_3.xdmf -O fine
```

</div>

<div class="w-full flex justify-center mt-2">
  <img src="./images/image_4.png" alt="Dirichlet Energy Minimization" class="w-full max-h-[180px] object-contain rounded shadow-sm" />
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

# 7-2. github 코드 실행 
generate_mesh_gmsh.py

::left::

<div class="text-sm mt-4">

#### 🔹 도메인 생성 (generate_disk)
* `gmsh.model.occ.addDisk`: 중심이 (0,0,0)이고 반지름이 1인 2D 원형 막(Disk)의 기하학적 형태를 생성합니다.<br><br>
#### 🔹 격자 세분화
* `mesh.refine()`: `refinement_level` 횟수만큼 격자를 반복해서 세분화합니다. 이를 통해 앞서 언급한 **메시 독립성(Mesh-independence) 검증**을 위한 다양한 해상도(coarse ~ fine)의 격자들을 준비합니다.<br><br>
#### 🔹 출력
* 맨 아래의 반복문(`for i in range(4)`)을 통해 4단계의 각기 다른 조밀도(레벨 0~3)를 가진 격자 파일들을 생성, 이후 최적화 솔버들이 불러와서 계산할 수 있도록 `.xdmf` 확장자로 데이터를 저장합니다.
</div>

::right::

<button @click="isExpanded = true" class="mt-8 px-4 py-2 bg-gray-800 text-white text-sm rounded shadow-md hover:bg-gray-700 transition-all flex items-center gap-2">
  <carbon:zoom-in /> 코드 전체화면으로 보기
</button>

<div :class="isExpanded ? 'fixed inset-4 z-50 bg-white shadow-2xl rounded-xl p-8 overflow-y-auto border border-gray-300' : 'overflow-y-auto max-h-[450px] shadow-lg rounded-md border border-gray-200/20 text-sm mt-8'">

  <div v-show="isExpanded" class="flex justify-between items-center mb-4 border-b pb-4">
    <div class="text-xl font-bold text-gray-800">generate_mesh_gmsh.py (전체 코드)</div>
    <button @click="isExpanded = false" class="px-4 py-1.5 bg-red-500 text-white rounded hover:bg-red-600 transition-all text-sm flex items-center gap-1">
      <carbon:close /> 닫기
    </button>
  </div>

  ```python {all} twoslash
  from pathlib import Path

  from mpi4py import MPI

  import dolfinx.io
  import gmsh
  import packaging.version

  __all__ = ["generate_disk"] 
  #generate_disk 외에 외부에서 접근 제한


  def generate_disk(filename: Path, res: float, order: int = 1, refinement_level: int = 1):
      """Generate a disk around the origin with radius 1 and resolution `res`.

      Args:
          filename: Name of the file to save the mesh to.
          res: Resolution of the mesh.
          order: Order of the mesh elements.
          refinement_level: Number of gmsh refinements
      """
      gmsh.initialize()
      if MPI.COMM_WORLD.rank == 0:
          membrane = gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
          gmsh.model.occ.synchronize()
          gdim = 2 #2D 문제
          gmsh.model.addPhysicalGroup(gdim, [membrane], 1)
          gmsh.option.setNumber("Mesh.CharacteristicLengthMin", res) 
          gmsh.option.setNumber("Mesh.CharacteristicLengthMax", res) #균일한 크기의 mesh 생성
          gmsh.model.mesh.generate(gdim)
          gmsh.model.mesh.setOrder(order)
          for _ in range(refinement_level):
              #refine level만큼 mesh 세분화(삼각형 1개 -> 4개)
              gmsh.model.mesh.refine()
              gmsh.model.mesh.setOrder(order)

      gmsh_model_rank = 0
      mesh_comm = MPI.COMM_WORLD
      model = dolfinx.io.gmsh.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)
      msh = model[0]
      gmsh.finalize()
      out_name = filename.with_stem(f"{filename.stem}_{refinement_level}").with_suffix(".xdmf")
      filename.parent.mkdir(exist_ok=True, parents=True)
      with dolfinx.io.XDMFFile(mesh_comm, out_name, "w") as xdmf:
          xdmf.write_mesh(msh)


  if __name__ == "__main__":
      for i in range(4):
          generate_disk(Path("meshes/disk.xdmf"), res=0.1, order=2, refinement_level=i)
  ```

</div>

---
layout: two-cols
layoutClass: gap-8
---

<script setup>
import { ref } from 'vue'
const isExpanded = ref(false)
</script>

# 7-3. github 코드 실행 
compare_all.py

::left::

<div class="text-sm mt-4">

#### 🔹 솔버 호출 (PG, IPOPT, GALAHAD, SNES)
* 각 코드의 솔버 함수를 import를 이용해 호출합니다.<br><br>
#### 🔹 출력 구성
* `Proximal Galerkin`: 1, 2차 라그랑주로 긱긱 출력
* `Galahad`: 단독 출력
* `IPOPT`: Hessian의 유무에 따른 두 경우 모두 출력
* `SNES`: 단독 출력
<div class="w-full flex flex-col items-center justify-center mt-6">
  <img src="./images/image_1.png" alt="Obstacle Cross Section" class="w-full max-h-[180px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-2">장애물 $\phi(x,y)$의 2D 단면 및 막의 접촉 형태</div>
</div>
</div>

::right::

<button @click="isExpanded = true" class="mt-8 px-4 py-2 bg-gray-800 text-white text-sm rounded shadow-md hover:bg-gray-700 transition-all flex items-center gap-2">
  <carbon:zoom-in /> 코드 전체화면으로 보기
</button>

<div :class="isExpanded ? 'fixed inset-4 z-50 bg-white shadow-2xl rounded-xl p-8 overflow-y-auto border border-gray-300' : 'overflow-y-auto max-h-[450px] shadow-lg rounded-md border border-gray-200/20 text-sm mt-8'">

  <div v-show="isExpanded" class="flex justify-between items-center mb-4 border-b pb-4">
    <div class="text-xl font-bold text-gray-800">generate_mesh_gmsh.py (전체 코드)</div>
    <button @click="isExpanded = false" class="px-4 py-1.5 bg-red-500 text-white rounded hover:bg-red-600 transition-all text-sm flex items-center gap-1">
      <carbon:close /> 닫기
    </button>
  </div>

  ```python {all} twoslash
  """
Solve obstacle problem with Proximal Galerkin, Galahad, and IPOPT and compare the results

Author: Jørgen S. Dokken
SPDX-License-Identifier: MIT
"""

import argparse
from pathlib import Path

import dolfinx
import numpy as np
from obstacle_ipopt_galahad import ObstacleProblem, setup_problem
from obstacle_pg import solve_problem
from obstacle_snes import snes_solve

from lvpp.optimization import galahad_solver, ipopt_solver

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Solve the obstacle problem on a unit square using Galahad.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--path",
        "-P",
        dest="infile",
        type=Path,
        default="./meshes/disk_3.xdmf",
        help="Path to infile",
    )
    parser.add_argument(
        "--results",
        "-O",
        dest="result_dir",
        type=Path,
        default=Path("results"),
        help="Path to results ",
    )
    max_iter = 500 # 최대 반복 횟수
    tol = 1e-4 # 수렴 허용 오차
    args = parser.parse_args()

    args.result_dir.mkdir(parents=True, exist_ok=True)

    # Set up problem matrices. initial guess and bounds
    problem = setup_problem(args.infile)
    S_, M_, f_, bounds_ = setup_problem(args.infile)
    # S_: 강성 행렬 (Stiffness matrix) — ∫∇u·∇v dx 에 해당
    # M_: 질량 행렬 (Mass matrix) — ∫u·v dx 에 해당
    # f_: 우변 하중 벡터 (FEM Function 형태)
    # bounds_: 장애물 조건 — [하한(lower), 상한(upper)] 형태의 Function 리스트

    bounds = tuple(b.x.array for b in bounds_)

    problem = ObstacleProblem(S_.copy(), M_.copy(), f_.x.array)
    # .copy()를 사용해 원본 행렬이 수정되지 않도록 보호
    V = f_.function_space
    mesh = V.mesh
    degree = mesh.geometry.cmap.degree
    V_out = dolfinx.fem.functionspace(mesh, ("Lagrange", degree))

    # Solve with Galahad
    x_g = dolfinx.fem.Function(V, name="galahad")
    x_g.x.array[:] = 0.0

    x_galahad, num_galahad_iterations = galahad_solver(
        problem,
        x_g.x.array.copy(),
        bounds,
        max_iter=max_iter,
        use_hessian=True,
        tol=tol,
    )
    x_g.x.array[:] = x_galahad
    x_g_out = dolfinx.fem.Function(V_out, name="ipopt")
    x_g_out.interpolate(x_g)
    with dolfinx.io.VTXWriter(
        V.mesh.comm, args.result_dir / f"{args.infile.stem}_galahad.bp", [x_g_out]
    ) as bp:
        bp.write(0.0)

    # Solve with llvp (first order)

    u_lvpp, max_it = solve_problem(
        args.infile,
        1, # 1차 라그랑주 사용
        maximum_number_of_outer_loop_iterations=max_iter,
        alpha_scheme="double_exponential",
        alpha_max=1e2, # alpha가 발산하지 않도록 최댓값 지정
        tol_exit=tol,
    )

    mesh = u_lvpp.function_space.mesh
    degree = mesh.geometry.cmap.degree
    V_out = dolfinx.fem.functionspace(mesh, ("Lagrange", degree))
    u_out = dolfinx.fem.Function(V_out, name="llvp")
    u_out.interpolate(u_lvpp.sub(0))
    with dolfinx.io.VTXWriter(
        mesh.comm, args.result_dir / f"{args.infile.stem}_llvp_first_order.bp", [u_out]
    ) as bp:
        bp.write(0.0)

    # Solve with llvp (second order)

    u_lvpp_2, max_it_2 = solve_problem(
        args.infile,
        2, # 2차 라그랑주 사용
        maximum_number_of_outer_loop_iterations=max_iter,
        alpha_scheme="double_exponential",
        alpha_max=1e2, # alpha가 발산하지 않도록 최댓값 지정
        tol_exit=tol,
    )
    u_out = u_lvpp_2.sub(0).collapse()
    with dolfinx.io.VTXWriter(
        u_out.function_space.mesh.comm,
        args.result_dir / f"{args.infile.stem}_llvp_second_order.bp",
        [u_out],
    ) as bp:
        bp.write(0.0)

    with dolfinx.io.VTXWriter(
        mesh.comm, args.result_dir / f"{args.infile.stem}_obstacle.bp", [bounds_[0]]
    ) as bp:
        bp.write(0.0)

    # Solve with IPOPT (With hessian)
    ipopt_iteration_count = {}
    for with_hessian in [True, False]:
        x_i = dolfinx.fem.Function(V, name="ipopt")
        x_i.x.array[:] = 0.0
        x_ipopt = ipopt_solver(
            problem,
            x_i.x.array.copy(),
            bounds,
            max_iter=max_iter,
            tol=1e-2 * tol,
            activate_hessian=with_hessian,
        )
        ipopt_iteration_count[with_hessian] = problem.total_iteration_count
        x_i.x.array[:] = x_ipopt

        # Output on geometry space

        x_i_out = dolfinx.fem.Function(V_out, name="ipopt")
        x_i_out.interpolate(x_i)
        with dolfinx.io.VTXWriter(
            mesh.comm,
            args.result_dir / f"{args.infile.stem}_ipopt_hessian_{with_hessian}.bp",
            [x_i_out],
        ) as bp:
            bp.write(0.0)

    # Solve with SNES
    u_snes, num_snes_iterations = snes_solve(
        args.infile,
        snes_options={
            "snes_type": "vinewtonssls",   # 변분 부등식용 Semi-smooth Newton
            "snes_monitor": None,
            "ksp_type": "preonly",          # 선형 솔버: 직접법(전처리만)
            "pc_type": "lu",               # 전처리기: LU 분해
            "pc_factor_mat_solver_type": "mumps",  # 병렬 직접 솔버 MUMPS 사용
        "snes_error_if_not_converged": True,   # 미수렴 시 예외 발생
            "ksp_error_if_not_converged": True,
        },
    )
    u_snes.name = "snes"
    with dolfinx.io.VTXWriter(
        mesh.comm,
        args.result_dir / f"{args.infile.stem}_snes.bp",
        [u_snes],
    ) as bp:
        bp.write(0.0)

    print(
        np.min(
            mesh.h(
                mesh.topology.dim, np.arange(mesh.topology.index_map(mesh.topology.dim).size_local)
            )
        )
    )
    print(f"{args.infile} Galahad iterations: {num_galahad_iterations}")
    print(f"{args.infile} llvp iterations: (P=1) {max_it}")
    print(f"{args.infile} llvp iterations: (P=2) {max_it_2}")
    print(f"{args.infile} Ipopt iterations: (With hessian) {ipopt_iteration_count[True]}")
    print(f"{args.infile} Ipopt iterations: (Without hessian {ipopt_iteration_count[False]}")
    print(f"{args.infile} SNES iterations: {num_snes_iterations}")
  ```

</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# GALAHAD — 내부점법 (Interior Point Method)
<div class="text-xl opacity-80 mb-6">로그 장벽 함수를 통한 제약 조건 내재화</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click>

### 🔹 핵심 아이디어
부등식 제약 $u \geq \phi$를 직접 다루는 대신, **로그 장벽 함수(Log Barrier)** 를 목적함수에 흡수시켜 제약 없는 문제로 변환합니다.

$$
J_\mu(u) = \frac{1}{2}u^TSu - f^TMu - \mu \sum_i \ln(u_i - \phi_i)
$$

<div class="mt-4 text-sm opacity-80">

- $\mu > 0$이 클수록 장벽이 강하게 작용하여 제약 경계 근처로의 접근을 억제합니다.
- $\mu \to 0$으로 점차 줄여가면 원래의 제약 문제 해에 수렴합니다.

</div>

### 🔹 수렴 흐름

$$
\underbrace{u_i \to \phi_i}_{\text{경계 접근}} \Rightarrow \ln(u_i - \phi_i) \to -\infty \Rightarrow J_\mu \to +\infty
$$

장벽 효과로 해가 **자연스럽게 실현 가능 영역(Feasible Region) 안에** 머뭅니다.

</div>
<div v-click>

### 🔹 알고리즘 단계

<div class="mt-2 text-sm opacity-80">

1. 초기 장벽 파라미터 $\mu_0 \gg 0$ 설정
2. 헤시안 $H = S + \mu \cdot \text{diag}((u_i - \phi_i)^{-2})$ 를 이용한 뉴턴 스텝으로 $J_\mu$ 최소화
3. $\mu_k \to 0$ 으로 감소시키며 반복
4. 수렴 조건 $\|Su - Mf\| < \text{tol}$ 달성 시 종료

</div><br>

### 🔹 헤시안 사용: 빠르게 최솟값에 도달

$$
u_{new}​=u_{old}​−H^{-1}∇J
$$

<div class="mt-2 text-sm opacity-80">

- `use_hessian=True` : 정확한 $H = S$ 사용 → **2차(Quadratic) 수렴**
- `use_hessian=False` : **L-BFGS** 근사 → 메모리 절약, 수렴 느림

</div>

</div>
</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# IPOPT — 내부점법 + 필터법
<div class="text-xl opacity-80 mb-6">KKT 조건의 직접 수치 해법</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click>

### 🔹 KKT 조건으로의 변환
제약 최적화 문제를 **KKT(Karush-Kuhn-Tucker) 조건** 으로 변환한 뒤, 이를 연립 비선형 방정식으로 풀어냅니다.

$$
\begin{cases}
Su - Mf - \lambda = 0 \quad \text{(정류 조건)}\\
(u - \phi) \cdot \lambda = \mu \quad \text{(상보성 조건)}\\
u - \phi \geq 0, \quad \lambda \geq 0
\end{cases}
$$

<div class="mt-2 text-sm opacity-80">

- $\lambda$: 제약에 대응하는 라그랑주 승수(반력)<br>
-> 막과 장애물이 떨어져 있다면 $\lambda=0$, 붙어 있다면 $\lambda \ge 0$
- $\mu \to 0$ 극한에서 정확한 상보성 조건 $(u-\phi)\lambda = 0$ 으로 수렴<br>
-> 막이 장애물과 닿아 있거나($u-\phi$ = 0) 떨어져 있거나($\lambda = 0$)
- 처음부터 막과 장애물이 닿아있는 부분을 모르므로 $\mu$를 양수로 한 후 천천히 줄여나감
</div>

</div>
<div v-click>

### 🔹 필터법 (Filter Method)
GALAHAD와의 핵심 차이점으로, 선탐색(Line Search)에 **필터법** 을 적용해 수렴 안정성을 높입니다.

<div class="mt-2 text-sm opacity-80">

- 목적함수 값과 제약 위반량을 동시에 추적하는 **파레토 필터** 유지
- 두 지표 중 하나라도 개선되는 스텝만 수락

</div><br>

### 🔹 헤시안 옵션 비교

| 옵션 | 방법 | 특징 |
|---|---|---|
| `activate_hessian=True` | 정확한 헤시안 $S$ | 빠른 수렴, 메모리 소모 |
| `activate_hessian=False` | **L-BFGS** 근사 | 느린 수렴, 메모리 절약 |

<div class="mt-2 text-sm opacity-80">

- 수렴 허용 오차를 다른 솔버 대비 **100배 더 엄격** 하게 설정
  (`tol = 1e-2 * tol`)

</div>

<div class="w-full flex justify-center mt-2">
  <img src="./images/image_2.png" alt="IPOPT Filter" class="w-full max-h-[120px] object-contain rounded shadow-sm border border-gray-200/50" />
</div>

</div>
</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# SNES — Semismooth Newton
<div class="text-xl opacity-80 mb-6">활성 집합을 이용한 변분 부등식의 직접 해법</div>

<div class="text-sm mt-4 mb-4">
<div grid="~ cols-2 gap-8" class="text-left">
<div v-click>

### 🔹 상보성 조건으로의 변환
장애물 문제의 최적성 조건(KKT)을 **상보성 조건** 으로 직접 표현합니다.

$$
\begin{cases}
u \geq \phi & \text{(실현 가능성)} \\
-\Delta u \geq 0 & \text{(쌍대 실현 가능성)} \\
(u - \phi)(-\Delta u) = 0 & \text{(상보성)}
\end{cases}
$$

<div class="mt-2 text-sm opacity-80">
- 막이 장애물에 **닿은 영역** : 반력 $-\Delta u > 0$이 발생
- 막이 장애물에서 **떨어진 영역** : 반력 = 0, $-\Delta u = 0$ (자유 평형)

</div><br>

### 🔹 활성 집합 (Active Set) 분류

$$
\mathcal{A} = \{i \mid u_i = \phi_i\}, \quad \mathcal{I} = \{i \mid u_i > \phi_i\}
$$

</div>
<div v-click>

### 🔹 vinewtonssls 알고리즘
PETSc의 `vinewtonssls`(Semismooth Line Search) 타입은 각 반복에서 활성 집합을 갱신하며 선형 시스템을 반복적으로 풀어냅니다.
-> 장애물과 접촉한 부분과 접촉하지 않은 부분을 갱신하며 영역을 둘로 나눠 계산.

<div class="mt-2 text-sm opacity-80">


</div><br>

### 🔹 경계 조건 설정 (코드)
```python
problem.solver.setVariableBounds(
    phi.x.petsc_vec,   # 하한: 장애물 φ
    u_max.x.petsc_vec  # 상한: +∞
)
```

<div class="mt-2 text-sm opacity-80">

- 별도의 장벽 파라미터 없이 **PDE 구조를 직접 보존**
- 선형 서브문제를 **LU 분해(MUMPS)** 로 정확하게 풀어 안정성 확보

</div>

</div>
</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 솔버 비교 요약
<div class="text-xl opacity-80 mb-6">GALAHAD · IPOPT · SNES의 접근 방식 차이</div>

<div grid="~ cols-2 gap-8" class="text-left text-sm mt-4">
<div v-click>

### 🔹 관점별 분류
```
장애물 문제 (부등식 제약 최적화)
│
├── GALAHAD / IPOPT ── 최적화 관점
│        │             "목적함수를 최소화하되
│        │              제약을 위반하지 말자"
│        └── 내부점법: 로그 장벽으로 제약을
│                      목적함수에 흡수
│
└── SNES ──────────── 방정식 관점
         │             "KKT 조건을 비선형
         │              방정식으로 직접 풀자"
         └── 활성집합법: 접촉 영역을
                        반복적으로 식별
```

</div>
<div v-click>

###

| | GALAHAD | IPOPT | SNES |
|---|---|---|---|
| **관점** | 최적화 | 최적화 | 비선형 방정식 |
| **핵심** | 로그 장벽 | 로그 장벽 + 필터 | 활성 집합 |
| **헤시안** | 선택 가능 | 선택 가능 | 항상 사용 |
| **선형 솔버** | 내장 | 내장 | MUMPS |
| **특징** | 범용 최적화 | 대규모에 강함 | PDE 특화 |

<div class="mt-2 text-sm opacity-80">

- **GALAHAD / SNES** : `tol = 1e-4`
- **IPOPT** : `tol = 1e-6` (100배 엄격)
- **공통** : `max_iter = 500`

</div>

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

# 7-4. github 코드 실행 
obstacle_pg.py

::left::

<div class="text-sm mt-4">

### 🔹 Mixed FunctionSpace
* 일반 솔버와 달리 물리 변수 $u$와 잠재 변수 $\psi$를 **동시에** 하나의 혼합 공간 $V = P_k \times P_k$ 에서 풉니다.<br><br>

### 🔹 비선형 잔차 $F$ 구성
* **$v$ 방정식**: $\alpha(\nabla u, \nabla v) + (\psi - \psi_k, v) = \alpha(f,v)$
* **$w$ 방정식** (latent): $u - e^\psi - \phi = 0$<br><br>

### 🔹 LVPP 외부 루프
* 매 반복마다 근위 파라미터 $\alpha_k$를 증가시키며 내부 SNES 서브문제를 반복 호출합니다.<br><br>

### 🔹 수렴 지표 기록
* 매 반복마다 에너지 · 상보성 · 실현가능성 · 뉴턴 스텝 수를 측정하여 `CSV`로 저장합니다.
</div>
::right::

<button @click="isExpanded = true" class="mt-8 px-4 py-2 bg-gray-800 text-white text-sm rounded shadow-md hover:bg-gray-700 transition-all flex items-center gap-2">
  <carbon:zoom-in /> 코드 전체화면으로 보기
</button>

<div :class="isExpanded ? 'fixed inset-4 z-50 bg-white shadow-2xl rounded-xl p-8 overflow-y-auto border border-gray-300' : 'overflow-y-auto max-h-[450px] shadow-lg rounded-md border border-gray-200/20 text-sm mt-8'">

  <div v-show="isExpanded" class="flex justify-between items-center mb-4 border-b pb-4">
    <div class="text-xl font-bold text-gray-800">generate_mesh_gmsh.py (전체 코드)</div>
    <button @click="isExpanded = false" class="px-4 py-1.5 bg-red-500 text-white rounded hover:bg-red-600 transition-all text-sm flex items-center gap-1">
      <carbon:close /> 닫기
    </button>
  </div>

  ```python {all} twoslash
  """

Obstacle problem based on experiment 4 in [2].

The FEniCSx code solve this problem is based on [3]:

SPXD License: MIT License

Original license file [../../licenses/LICENSE.surowiec](../../licenses/LICENSE.surowiec)
is included in the repository.

[1] Keith, B. and Surowiec, T.M., Proximal Galerkin: A Structure-Preserving Finite Element Method
for Pointwise Bound Constraints. Found Comput Math (2024). https://doi.org/10.1007/s10208-024-09681-8
[2] Keith, B., Surowiec, T. M., & Dokken, J. S. (2023). Examples for the Proximal Galerkin Method
    (Version 0.1.0) [Computer software]. https://github.com/thomas-surowiec/proximal-galerkin-examples
"""

import argparse
from pathlib import Path

from mpi4py import MPI

import basix
import numpy as np
import pandas as pd
import ufl
from dolfinx import default_scalar_type, fem, io, mesh
from dolfinx.fem.petsc import NonlinearProblem
from ufl import Measure, conditional, exp, grad, inner, lt


def rank_print(string: str, comm: MPI.Comm, rank: int = 0):
    """Helper function to print on a single rank

    :param string: String to print
    :param comm: The MPI communicator
    :param rank: Rank to print on, defaults to 0
    """
    if comm.rank == rank:
        print(string)


def allreduce_scalar(form: fem.Form, op: MPI.Op = MPI.SUM) -> np.floating:
    """Assemble a scalar form over all processes and perform a global reduction

    :param form: Scalar form
    :param op: MPI reduction operation
    """
    comm = form.mesh.comm
    return comm.allreduce(fem.assemble_scalar(form), op=op)


def solve_problem(
    filename: Path,
    polynomial_order: int,
    maximum_number_of_outer_loop_iterations: int,
    alpha_scheme: str,
    alpha_max: float,
    tol_exit: float,
):
    """ """

    # Create mesh
    with io.XDMFFile(MPI.COMM_WORLD, filename, "r") as xdmf:
        msh = xdmf.read_mesh(name="mesh")

    # Define FE subspaces
    P = basix.ufl.element("Lagrange", msh.basix_cell(), polynomial_order)
    mixed_element = basix.ufl.mixed_element([P, P])
    V = fem.functionspace(msh, mixed_element)

    # Define functions and parameters
    alpha = fem.Constant(msh, default_scalar_type(1))
    f = fem.Constant(msh, 0.0)
    # Define BCs
    msh.topology.create_connectivity(msh.topology.dim - 1, msh.topology.dim)
    facets = mesh.exterior_facet_indices(msh.topology)
    V0, _ = V.sub(0).collapse()
    dofs = fem.locate_dofs_topological((V.sub(0), V0), entity_dim=1, entities=facets)

    u_bc = fem.Function(V0)
    u_bc.x.array[:] = 0.0
    bcs = fem.dirichletbc(value=u_bc, dofs=dofs, V=V.sub(0))

    # Define solution variables
    sol = fem.Function(V)
    sol_k = fem.Function(V)

    u, psi = ufl.split(sol) # 두 성분을 각각 분리하여 계산에 이용
    u_k, psi_k = ufl.split(sol_k)

    def phi_set(x): # 장애물 생성 함수
        r = np.sqrt(x[0] ** 2 + x[1] ** 2)
        r0 = 0.5 # 반지름 = 0.5
        beta = 0.9
        b = r0 * beta
        tmp = np.sqrt(r0**2 - b**2)
        B = tmp + b * b / tmp
        C = -b / tmp
        cond_true = B + r * C
        cond_false = np.sqrt(r0**2 - r**2)
        true_indices = np.flatnonzero(r > b)
        cond_false[true_indices] = cond_true[true_indices]
        return cond_false

    quadrature_degree = 6
    Qe = basix.ufl.quadrature_element(msh.topology.cell_name(), degree=quadrature_degree)
    Vq = fem.functionspace(msh, Qe)
    # Lower bound for the obstacle
    phi = fem.Function(Vq, name="phi")
    phi.interpolate(phi_set)

    # Define non-linear residual
    (v, w) = ufl.TestFunctions(V)
    dx = Measure("dx", domain=msh, metadata={"quadrature_degree": quadrature_degree})
    F = (
        alpha * inner(grad(u), grad(v)) * dx   # ① α(∇u, ∇v)
        + psi * v * dx                          # ② (ψ, v)
        + u * w * dx                            # ③ (u, w)
        - exp(psi) * w * dx                     # ④ -(e^ψ, w)
        - phi * w * dx                          # ⑤ -(φ, w)
        - alpha * f * v * dx                    # ⑥ -α(f, v)
        - psi_k * v * dx                        # ⑦ -(ψ_k, v)
    )
    J = ufl.derivative(F, sol)

    # Setup non-linear problem
    petsc_options = {
        "ksp_type": "preonly", # 선형 솔버: 직접법
        "pc_type": "lu", # 전처리: LU 분해
        "pc_factor_mat_solver_type": "mumps", # 병렬 직접 솔버 MUMPS
        "ksp_error_if_not_converged": True,
        "ksp_monitor": None,
        "snes_monitor": None,
        "snes_error_if_not_converged": True,
        "snes_linesearch_type": "none", # 선탐색 없음 (순수 뉴턴)
        "snes_rtol": 1e-6, # 상대 수렴 허용 오차
        "snes_max_it": 100, # 최대 뉴턴 반복 횟수     
    }
    problem = NonlinearProblem(
        F, u=sol, bcs=[bcs], J=J, petsc_options=petsc_options, petsc_options_prefix="obstacle_"
    )

    # observables
    energy_form = fem.form(0.5 * inner(grad(u), grad(u)) * dx - f * u * dx) # 에너지 식.
    complementarity_form = fem.form((psi_k - psi) / alpha * u * dx) # 반복 계산 간 차이
    feasibility_form = fem.form(conditional(lt(u, 0), -u, fem.Constant(msh, 0.0)) * dx) # 음수 영역의 크기를 측정(제약 위반인지 확인)
    dual_feasibility_form = fem.form(
        conditional(lt(psi_k, psi), (psi - psi_k) / alpha, fem.Constant(msh, 0.0)) * dx
    ) # ψ_k - ψ >= 0인지 확인
    H1increment_form = fem.form(inner(grad(u - u_k), grad(u - u_k)) * dx + (u - u_k) ** 2 * dx)
    # u의 변화량 확인 후 이 값이 tol_exit보다 작으면 계산 종료
    L2increment_form = fem.form((exp(psi) - exp(psi_k)) ** 2 * dx)
    # ψ의 변화량을 e^ψ, 즉 u의 영역에서 측정

    # Proximal point outer loop
    n = 0
    increment_k = 0.0
    sol.x.array[:] = 0.0
    sol_k.x.array[:] = sol.x.array[:]
    alpha_k = 1 #근위 파라미터, 반복마다 업데이트되는 스칼라 상수
    step_size_rule = alpha_scheme
    C = 1.0
    r = 1.5
    q = 1.5

    energies = []
    complementarities = []
    feasibilities = []
    dual_feasibilities = []
    Newton_steps = []
    step_sizes = []
    primal_increments = []
    latent_increments = []
    for k in range(maximum_number_of_outer_loop_iterations):
        # Update step size
        if step_size_rule == "constant":
            alpha.value = C # α = 1 고정
        elif step_size_rule == "double_exponential":
            try:
                alpha.value = max(C * r ** (q**k) - alpha_k, C) # α ~ 1.5^(1.5^k)
            except OverflowError:
                pass
            alpha_k = alpha.value
            alpha.value = min(alpha.value, alpha_max) # 상한 적용
        else:
            step_size_rule == "geometric"
            alpha.value = C * r**k # α = 1.5^k (등비수열)
        rank_print(f"OUTER LOOP {k + 1} alpha: {alpha.value}", msh.comm)

        # Solve problem
        problem.solve()
        converged_reason = problem.solver.getConvergedReason()
        n = problem.solver.getIterationNumber()
        rank_print(f"Newton steps: {n}   Converged: {converged_reason}", msh.comm)

        # Check outer loop convergence
        energy = allreduce_scalar(energy_form)
        complementarity = np.abs(allreduce_scalar(complementarity_form))
        feasibility = allreduce_scalar(feasibility_form)
        dual_feasibility = allreduce_scalar(dual_feasibility_form)
        increment = np.sqrt(allreduce_scalar(H1increment_form))
        latent_increment = np.sqrt(allreduce_scalar(L2increment_form))

        tol_pp = increment

        if increment_k > 0.0:
            rank_print(
                f"Increment size: {increment}" + f"   Ratio: {increment / increment_k}", msh.comm
            )
        else:
            rank_print(f"Increment size: {increment}", msh.comm)
        rank_print("", msh.comm)

        energies.append(energy)
        complementarities.append(complementarity)
        feasibilities.append(feasibility)
        dual_feasibilities.append(dual_feasibility)
        Newton_steps.append(n)
        step_sizes.append(np.copy(alpha.value))
        primal_increments.append(increment)
        latent_increments.append(latent_increment)

        if tol_pp < tol_exit:
            break

        # Update sol_k with sol_new
        sol_k.x.array[:] = sol.x.array[:]
        increment_k = increment

    # # Save data
    cwd = Path.cwd()
    output_dir = cwd / "output"
    output_dir.mkdir(exist_ok=True)

    # Create output space for bubble function
    V_primal, primal_to_mixed = V.sub(0).collapse()

    num_primal_dofs = V_primal.dofmap.index_map.size_global

    phi_out_space = fem.functionspace(msh, basix.ufl.element("Lagrange", msh.basix_cell(), 6))
    phi_out = fem.Function(phi_out_space, name="phi")
    phi_out.interpolate(phi_set)
    with io.VTXWriter(msh.comm, output_dir / "phi.bp", [phi_out]) as bp:
        bp.write(0.0)
    if MPI.COMM_WORLD.rank == 0:
        df = pd.DataFrame()
        df["Energy"] = energies
        df["Complementarity"] = complementarities
        df["Feasibility"] = feasibilities
        df["Dual Feasibility"] = dual_feasibilities
        df["Newton steps"] = Newton_steps
        df["Step sizes"] = step_sizes
        df["Primal increments"] = primal_increments
        df["Latent increments"] = latent_increments
        df["Polynomial order"] = np.full(k + 1, polynomial_order)
        df["dofs"] = np.full(k + 1, num_primal_dofs)
        df["Step size rule"] = [step_size_rule] * (k + 1)
        filename = f"./example_polyorder{polynomial_order}_{num_primal_dofs}.csv"
        print(f"Saving data to: {str(output_dir / filename)}")
        df.to_csv(output_dir / filename, index=False)
        rank_print(df, msh.comm)

    if k == maximum_number_of_outer_loop_iterations - 1:
        rank_print("Maximum number of outer loop iterations reached", msh.comm)
    return sol, sum(Newton_steps)


# -------------------------------------------------------
if __name__ == "__main__":
    desc = "Run examples from paper"
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--file-path",
        "-f",
        dest="filename",
        type=Path,
        required=True,
        help="Name of input file",
    )
    parser.add_argument(
        "--polynomial_order",
        "-p",
        dest="polynomial_order",
        type=int,
        default=1,
        choices=[1, 2],
        help="Polynomial order of primal space",
    )
    parser.add_argument(
        "--alpha-scheme",
        dest="alpha_scheme",
        type=str,
        default="constant",
        choices=["constant", "double_exponential", "geometric"],
        help="Step size rule",
    )
    parser.add_argument(
        "--max-iter",
        "-i",
        dest="maximum_number_of_outer_loop_iterations",
        type=int,
        default=100,
        help="Maximum number of outer loop iterations",
    )
    parser.add_argument(
        "--alpha-max",
        "-a",
        dest="alpha_max",
        type=float,
        default=1e5,
        help="Maximum alpha",
    )
    parser.add_argument(
        "--tol",
        "-t",
        dest="tol_exit",
        type=float,
        default=1e-6,
        help="Tolerance for exiting Newton iteration",
    )
    args = parser.parse_args()
    solve_problem(
        args.filename,
        args.polynomial_order,
        args.maximum_number_of_outer_loop_iterations,
        args.alpha_scheme,
        args.alpha_max,
        args.tol_exit,
    )
  ```

</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 8-1. 실행 결과 분석

<div class="text-xl opacity-80 mb-8">
  각 솔버의 수치 해석 결과 비교
</div>

<div grid="~ cols-2 gap-8" class="text-left">

<div v-click>

### 🔹 GALAHAD
<div class="w-full flex flex-col items-center">
  <img src="./images/image_2.png" alt="GALAHAD 결과" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">GALAHAD 솔버 실행 결과</div>
</div>

</div>

<div v-click>

### 🔹 SNES
<div class="w-full flex flex-col items-center">
  <img src="./images/image_2.png" alt="IPOPT 결과" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">SNES 솔버 실행 결과</div>
</div>

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 8-2. 실행 결과 분석

<div class="text-xl opacity-80 mb-8">
  각 솔버의 수치 해석 결과 비교
</div>

<div grid="~ cols-2 gap-8" class="text-left">

<div v-click>

### 🔹 IPOPT(with hessian)
<div class="w-full flex flex-col items-center">
  <img src="./images/image_2.png" alt="GALAHAD 결과" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">IPOPT(with hessian) 솔버 실행 결과</div>
</div>

</div>

<div v-click>

### 🔹 IPOPT(without hessian)
<div class="w-full flex flex-col items-center">
  <img src="./images/image_2.png" alt="IPOPT 결과" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">IPOPT(without hessian) 솔버 실행 결과</div>
</div>

</div>
</div>
---
theme: seriph
class: text-center
highlighter: shiki
---

# 8-3. 실행 결과 분석

<div class="text-xl opacity-80 mb-8">
  각 솔버의 수치 해석 결과 비교
</div>

<div grid="~ cols-2 gap-8" class="text-left">

<div v-click>

### 🔹 Proximal Galerkin(1차)
<div class="w-full flex flex-col items-center">
  <img src="./images/image_2.png" alt="GALAHAD 결과" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">Proximal Galerkin(1차) 솔버 실행 결과</div>
</div>

</div>

<div v-click>

### 🔹 Proximal Galerkin(2차)
<div class="w-full flex flex-col items-center">
  <img src="./images/image_2.png" alt="IPOPT 결과" class="w-full max-h-[160px] object-contain rounded shadow-sm border border-gray-200/50" />
  <div class="text-xs opacity-70 mt-1">Proximal Galerkin(2차) 솔버 실행 결과</div>
</div>

</div>
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 9. 향후 발전 방향 (Develop Direction)
<div class="text-xl opacity-80 mb-6">멤브레인을 이용한 고체 슬러리/액체 분리 모델링</div>

<div class="text-sm">
<div grid="~ cols-2 gap-4" class="mt-4 mb-4">
<div>

### 🔹 공학적 응용 프레임워크 확장
* 본 알고리즘이 적용된 고체 접촉(Membrane-Obstacle) 모델을 응용하여 **화공생명공학의 분리 공정 모델링**으로 발전시킬 수 있습니다.
* 여과막(Membrane)을 통과하는 액체와 통과하지 못하는 고체 슬러리(Slurry) 혼합물 분리 시스템에 적용합니다.

</div>
<div>

### 🔹 점별 값 제약(Value Constraint)의 활용
* 고체 입자가 여과막을 뚫지 못하고 쌓이면서 Cake layer를 형성하는 물리적 현상을 $u \ge \phi$ (입자 위치 $\ge$ 여과막 표면) 형태의 제약으로 모사합니다.
* 이를 FEniCSx 기반 LVPP 솔버로 계산하면, 격자 의존성 없이 압력 분포 및 막의 미세 변형을 대규모로 병렬 처리할 수 있을 것으로 기대됩니다.

</div>
</div>
</div>

<div class="w-full flex justify-center mt-2">
  <img src="./images/image_1.png" alt="Slurry Separation Model" class="w-full max-h-[220px] object-contain rounded shadow-sm" />
</div>

---
layout: center
class: text-center
---

# Thank You

<div class="mt-8 text-xl">
  Q & A
</div>
