---
theme: seriph
title: 4월 7일 팀 FEniCSx 발표
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

# 4월 7일 팀 FEniCSx 발표

<div class="text-xl opacity-80 mb-8 mt-4">
  논문 예제의 Developing 과정 및 논문 준비 현황
</div>

</div>

<div class="absolute bottom-6 right-6 text-xl z-10" color='black'>
  Team FEniCSx: 박기성 / 윤현준 / 이민용 / 박찬서
</div>

---
theme: seriph
class: text-center
highlighter: shiki
---

# 논문 준비 현황

<div class="text-left mt-12 text-lg space-y-5 px-10">

<div v-click>

1. **수식 정리 완료**
나비에-스토크스(Navier-Stokes) 및 막(Membrane) 거동에 대한 논문용 수식 정리를 마쳤습니다.

</div>
<div v-click>

2. **Develop 아이디어 구상**
논문을 발전시킬 후보 주제 3가지를 도출했습니다. (어떤 것을 최종 적용할지는 아직 미정입니다.)

</div>
<div v-click>

3. **시뮬레이션 진행 상황**
FEniCSx를 이용해 지속적으로 시뮬레이션을 시도 중이나, 아직 뚜렷하게 성공한 결과는 나오지 않았습니다.

</div>
<div v-click>

4. **수치적 한계점 파악**
시뮬레이션을 진행하면서 극단적으로 많은 타임스텝 수(dt)가 요구되는 등 컴퓨팅 측면의 한계점을 분석하고 있습니다.

</div>
<div v-click>

5. **논문 목차 구성**
서론-본론-결론으로 이어지는 전반적인 논문의 흐름과 뼈대(목차) 구성을 완료했습니다.

</div>
</div>

---
theme: seriph
class: text-center
---

# 논문 목차 (Draft)
<div class="text-xl opacity-80 mb-8">4단 논리 흐름 구성</div>

<div class="grid grid-cols-2 gap-8 text-left">
  
<div>

### 1. 서론 (Introduction)
* 화공 공정 설계모사 내 국소적 유동 모사의 필요성
* LVPP 알고리즘 기반 검증 모델 소개 (예제 1, 5)

<br>

### 2. 이론 및 해석 모델 (Theoretical Background)
* **지배 방정식 및 물리적 구조별 접근:**
  * 유동 영역: Navier-Stokes 방정식
  * 다공성 막(Membrane) 영역: Darcy 법칙
  * 메쉬 변형(Mesh update): ALE 기법 및 Poisson 방정식

</div>

<div>

### 3. 결과 및 고찰 (Results & Discussion)
* **시뮬레이션 결과 분석:** 막의 변형 거동 및 영역별 유동 프로파일 검증
* **연구의 확장성:** 제안 모델의 산업적 적용 가치 및 필요성

<br>

### 4. 결론 (Conclusion)
* 기초 모델의 산업적 확대 적용 시 기대 효과
* 한계점 및 향후 연구 방향

</div>

</div>

---
layout: two-cols
layoutClass: gap-8
---

# 현 시뮬레이션에서의 어려움 및 현황
<div class="text-xl opacity-80 mb-8">메시 구성의 어려움</div>

::left::

<div class="text-[14px] mt-4 space-y-4 pr-4">

### 🔹 주요 이슈
- 메시 크기가 커서 장애물과 접촉하는 부분을 알기 어렵습니다.(1번 develop)
- 관-막의 시뮬레이션에서 LVPP를 적용할 때 장애물 메시 구성에 어려움이 있습니다.(2번 develop)

### 🔹 어떻게 해결할 것인지.
- 막의 장력 등 파라미터 조정해서 시각화를 용이하게 하기 위해 시도중입니다.(1번 develop)
- 장애물 메시를 원관 메시와 독립적으로 구성해 lvpp 적용할 생각입니다.(2번 develop)

</div>

::right::

<div class="overflow-y-auto max-h-[400px] shadow-lg rounded-md border border-gray-200/20 text-sm mt-4">

```python {all} twoslash
FSI: ALE Navier-Stokes + Darcy + Poisson Membrane + Sphere Obstacle (LVPP)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
전체 연성 구조:
  ① NS IPCS (유체 속도·압력) + Brinkman 항 (장애물 내부 유동 차단)
  ② 구 운동 방정식 (Stokes 항력 + 접촉 법선력)
  ③ chi_s 갱신 (구 이동 후 Brinkman 영역 재계산)
  ④ φ 갱신 (구 위치 기반 막 장애물 함수)
  ⑤ pressure_jump_dynamic (막 전후 압력차 계산 → LVPP 우변)
  ⑥ LVPP 막 (w ≥ φ 보장, 혼합 공간 (w,ψ))
  ⑦ ALE 조화 확장 (막 변위 → 유체 메시 변형)
```
</div>

---
theme: seriph
class: text-center
---

# 앞으로 할 것 (Future Works)

<div class="text-left mt-16 text-xl space-y-10 px-16">

<div v-click>

### 1. 시뮬레이션 시도 후 성공 시 연구실 컴퓨터 활용

<br>

- lvpp가 작동하는 것을 확인 시 메시를 조밀하게 하여 시각화에 용이하도록 만들어 논문에 사용하도록 연구실 컴퓨터 사용할 예정입니다.

</div>

<div v-click>

### 2. 관련 문헌 집중 조사

<br>

현재 구축 중인 논리 구조의 타당성을 확보하기 위해 관련 논문들을 추가로 읽고 정리할 계획입니다.
* ALE (Arbitrary Lagrangian-Eulerian) 기법을 적용한 상호작용 논문
* 멤브레인 필터(Membrane Filter) 구조 확인

</div>
</div>
---
layout: center
class: text-center
---

# Thank You

<div class="mt-8 text-xl">
  Q & A
</div>