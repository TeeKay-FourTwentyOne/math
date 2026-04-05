# Iteration 11: Cross-Platform Noise Calibration

**Direction**: Calibrate which row of our divergence table applies to real hardware differences, using published empirical data on FP noise in GPU computations.

---

## Literature Findings

### 1. GEMM-Level Noise (same GPU, batch-dependent)
**Source**: "On the Structure of Floating-Point Noise in Batch-Invariant GPU Matrix Multiplication" (arXiv 2511.00025)
- **BF16**: sigma = 1.17 x 10^-3 per element
- **FP16**: sigma = 5.32 x 10^-4 per element
- Noise is **structured and correlated** (47% of FP16 variance is off-diagonal)
- A **prediction flip** was demonstrated at absolute deviation 2 x 10^-4 with logit margin 0.02

### 2. Cross-GPU GEMM Noise (different NVIDIA hardware)
**Source**: Ingonyama engineering blog ("Solving Reproducibility Challenges")
- GEMM outputs differed by ~**1e-4** between L4, RTX 3090, and RTX 4080
- Quote: "the results differed with an error margin of 1e-4"
- Note: these are all FP32 operations on CUDA GPUs of different architectures

### 3. PyTorch Operation-Level Variability (FP32)
**Source**: "Impacts of floating-point non-associativity on reproducibility" (SC'24, arXiv 2408.05148)
- scatter_reduce: variability 0 to 3.35 x 10^-6
- index_add: variability 0 to 5.03 x 10^-6
- These are individual FP32 operations on H100

### 4. Same-GPU Batch Nondeterminism
**Source**: Thinking Machines Lab blog ("Defeating Nondeterminism in LLM Inference")
- 1000 completions at temperature 0: **80 unique outputs**
- First divergence at token **103**
- Precision not specified (likely FP16/BF16 based on production deployment)

## Accumulation Model

Cross-GPU logit noise accumulates through the transformer:

**Per-GEMM noise estimate** (FP32, cross-GPU): ~1e-4 per element (Ingonyama)

**Accumulation through 24 layers**:
- ~4 GEMMs per layer (Q, K, V projections + MLP up/down)
- 96 total GEMMs per forward pass
- Independent accumulation (sqrt model): sqrt(96) x 1e-4 ~ **1e-3**
- Worst-case accumulation (linear model): 96 x 1e-4 ~ **1e-2**
- Realistic estimate: between sqrt and linear, call it **1e-3 to 5e-3**

**Additional factor**: autoregressive generation compounds errors across 256 steps. Each step's output feeds the next step, so errors can grow. However, the KV cache means prior layer computations are not recomputed — only the current token's forward pass introduces new errors.

## Calibrated Predictions

| Platform Comparison | Precision | Estimated logit epsilon | Our predicted divergence rate |
|-------------------|-----------|----------------------|------------------------------|
| Same GPU, deterministic | FP32 | 0 (exact) | **0%** (verified) |
| Same GPU, non-deterministic | FP32 | ~1e-6 to 1e-5 | **~0-0.05%** |
| **Different NVIDIA GPUs** | **FP32** | **~1e-4 to 1e-3** | **~0.3-4.4%** |
| **NVIDIA vs Apple Silicon** | **FP32** | **~5e-4 to 5e-3** | **~2-20%** |
| Any cross-platform | BF16 | ~1e-3 to 1e-2 | **~4-32%** |
| Any cross-platform | FP16 | ~5e-4 to 5e-3 | **~2-20%** |

### Interpreting for Our Experiment

For Experiment 1A (same model, different platforms, FP32, greedy decoding):

**H_nvidia (RTX 4070) vs H_cloud (A100)**: Both NVIDIA CUDA, but different architectures
- Expected epsilon: ~1e-4 to 1e-3
- Predicted: **0.3-4.4% of 400 prompts diverge** (1-18 prompts)
- Experimental prompts: 9x more likely to diverge than control

**H_nvidia (RTX 4070) vs H_mac (Apple Silicon)**: Different GPU vendors, different ISAs
- Expected epsilon: ~5e-4 to 5e-3 (larger due to fundamental architecture differences)
- Predicted: **2-20% of prompts diverge** (8-80 prompts)
- This is the most architecturally interesting comparison

**For Experiment 1C (BF16/FP16)**:
- Expected epsilon: 10-100x larger than FP32
- Predicted: **20-75% of prompts diverge** — massive divergence expected
- This should be the most dramatic result in the paper

## Key Insight: FP32 May Be "Too Clean"

The literature suggests that FP32 cross-GPU noise (~1e-4 per GEMM) may only flip tokens at the very tightest margins. Our min margin across all prompts is 0.0004 — right at the threshold. This means:

- FP32 cross-platform divergence will be **sparse but real**: only the most fragile positions will flip
- BF16/FP16 will show **dramatic divergence**: many more positions have margins < 1e-3
- The precision comparison (Experiment 1C) may be more revealing than the platform comparison (1A)

## Caveats

1. The 1e-4 cross-GPU GEMM noise is from comparing **different NVIDIA architectures** (Ampere vs Turing vs Ada). NVIDIA vs Apple Silicon could be worse.
2. Error accumulation through 24 layers depends on whether errors are correlated or independent. The 2511.00025 paper showed errors are **highly structured and correlated** (not independent Gaussian), which could mean faster-than-sqrt accumulation.
3. Our experiment uses torch.use_deterministic_algorithms(True), which eliminates SOME sources of within-GPU noise but NOT cross-GPU noise (different FP implementations produce deterministically different results).
4. Autoregressive compounding means later tokens accumulate more error — consistent with our finding that late positions show more divergence.

## Sources
- [On the Structure of Floating-Point Noise in Batch-Invariant GPU Matrix Multiplication](https://arxiv.org/html/2511.00025)
- [Solving Reproducibility Challenges in Deep Learning and LLMs](https://www.ingonyama.com/post/solving-reproducibility-challenges-in-deep-learning-and-llms-our-journey)
- [Impacts of floating-point non-associativity on reproducibility (SC'24)](https://arxiv.org/html/2408.05148v3)
- [Defeating Nondeterminism in LLM Inference](https://thinkingmachines.ai/blog/defeating-nondeterminism-in-llm-inference/)

## Scripts
- `research_loop/iter11_noise_calibration.md` (this report — no code needed, literature synthesis)
