# Experimental Protocol: Hardware-Induced Divergence in LLM Training and Inference

## Overview

This document describes a set of experiments to empirically measure how different hardware produces different outputs from neural language models, isolating two independent sources of divergence:

- **Scenario 1 (Inference divergence):** Same weights, same prompt, different hardware. Measures how floating-point implementation differences during inference cause output divergence.
- **Scenario 2 (Training divergence):** Same training recipe but different training hardware, producing models with slightly different weights. Measures how those weight differences compound with inference divergence.

The experiments use a shared prompt set, shared measurement framework, and a 2x2 factorial design (same/different training hardware × same/different inference hardware) to cleanly separate variance sources.

---

## Hardware Platforms

Three platforms are used throughout:

| Label | Hardware | Notes |
|-------|----------|-------|
| **H_nvidia** | Windows desktop, NVIDIA consumer GPU | CUDA backend. Record exact GPU model and driver version. |
| **H_mac** | MacBook, Apple Silicon | MPS backend (or CPU fallback). Record chip model (M1/M2/M3). |
| **H_cloud** | Cloud instance, NVIDIA datacenter GPU | Aim for A100 or L40S. RunPod or Lambda recommended. ~$1-2/hr. |

Before any experiment, record and save a hardware manifest for each platform: GPU model, driver version, PyTorch version, Python version, OS version, CUDA/MPS version. This goes into a `hardware_manifest.json` for each platform.

---

## Model Selection

**Primary model:** Pythia-410M (EleutherAI/pythia-410m-deduped)

Rationale:
- Designed for reproducibility research; training data and procedure fully documented
- 410M parameters, 24 transformer layers — small enough to record full activations, large enough to generate coherent text
- Runs on all three platforms without quantization

**Fallback model (if 410M is too slow or large for activation recording):** Pythia-160M or GPT-2 Small (124M, 12 layers). Smaller means faster activation dumps and simpler divergence profiles.

**Precision configurations to test:** FP32 (baseline, most deterministic), BF16 (if supported on all platforms), FP16 (expected to show more divergence than FP32).

---

## Experiment 0: Prompt Set Construction (Run First)

### Goal
Select prompts that maximize the probability of detecting hardware-induced divergence. Prompts where the model is near a token-selection decision boundary are most informative.

### Procedure

1. Assemble a candidate pool of ~2,000 prompts across five categories:
   - **Factual recall (400):** "What is the capital of France?", "Name three elements heavier than iron." These are expected to be robust (high-confidence answers). They serve as negative controls.
   - **Open-ended generation (400):** "Write a paragraph about the feeling of rain on a tin roof.", "Describe a city that doesn't exist." Expected to be fragile.
   - **Reasoning/math (400):** "What is 247 + 389?", "If a train leaves at 3pm traveling 60mph..." Mixed fragility.
   - **Continuation (400):** Provide 2-3 sentences of text, ask the model to continue. Use passages from public domain sources. Expected fragile.
   - **Conversational (400):** "Hey, how's it going?", "What do you think about pineapple on pizza?" Expected moderate fragility.

2. Run all 2,000 prompts on ONE platform (H_nvidia recommended) with greedy decoding, temperature=0, generating 256 tokens per prompt.

3. At every token generation step, log the **top-5 logits** (or at minimum, the logit values for the rank-1 and rank-2 tokens).

4. For each prompt, compute a **fragility score:**
   ```
   fragility = (number of positions where |logit_rank1 - logit_rank2| < threshold) / total_positions
   ```
   Start with threshold = 1.0 in logit space. Adjust if too many or too few prompts qualify.

5. Also compute **minimum margin** (smallest gap between top-2 logits across the sequence) and **mean margin** per prompt.

6. Rank prompts by fragility score. Select:
   - **Top 200 most fragile** prompts (the experimental set — maximizes signal)
   - **200 least fragile** prompts (the control set — expected to show no divergence)
   - Ensure both sets have roughly equal representation across the 5 categories

7. Save the selected prompt sets, their fragility scores, and full logit traces as the fixed experimental inputs for all subsequent experiments.

### Output Artifacts
- `prompts_experimental.jsonl` — 200 high-fragility prompts with metadata
- `prompts_control.jsonl` — 200 low-fragility prompts with metadata
- `fragility_scores.csv` — all 2,000 prompts with scores
- `logit_traces/` — directory of per-prompt logit logs from the scoring run

---

## Experiment 1: Inference Divergence (Scenario 1)

### Goal
Measure output divergence when identical weights process identical prompts on different hardware.

### Setup
- Model: Pythia-410M, identical checkpoint downloaded to all three platforms
- Verify checkpoints are bitwise identical (SHA-256 hash of all weight files must match)
- Prompts: The 400-prompt set from Experiment 0 (200 experimental + 200 control)

### Procedure

**Phase 1A — Token-level divergence measurement:**

For each platform (H_nvidia, H_mac, H_cloud):
1. Load model in FP32
2. Run all 400 prompts, greedy decoding, temperature=0, max_new_tokens=256
3. Record full token ID sequences
4. Repeat the entire run a second time on the same hardware (to measure within-platform nondeterminism baseline)

This produces 6 output sets (3 platforms × 2 runs each). Compare:
- Within-platform: Run 1 vs Run 2 on same hardware (baseline nondeterminism)
- Cross-platform: Each pairwise comparison (nvidia vs mac, nvidia vs cloud, mac vs cloud)

**Phase 1B — Activation-level divergence profiling:**

For a subset of 50 prompts (25 fragile, 25 control):
1. On each platform, run inference with forward hooks that record the output tensor of every transformer layer
2. Record at full precision (FP32 regardless of model precision)
3. For each prompt, for each layer, compute:
   - Element-wise absolute difference between platforms
   - L2 norm of the difference
   - Max absolute difference
   - Relative difference (normalized by activation magnitude)

This produces a **divergence profile** per prompt: a curve showing how error magnitude evolves layer by layer.

**Phase 1C — Precision sensitivity:**

Repeat Phase 1A at BF16 and FP16 (if supported on all platforms). Compare divergence rates across precisions.

### Measurements to Record

For each pair of output sequences (per prompt):
- **Exact match:** Boolean — are all 256 tokens identical?
- **First divergence point:** Token index where outputs first differ (or null if exact match)
- **Token agreement rate:** Fraction of 256 positions with identical tokens
- **Margin at divergence point:** What was the logit gap at the position where divergence occurred? (from Experiment 0 logit traces)
- **Levenshtein distance** of detokenized text
- **Cosine similarity** of sentence embeddings (use a small embedding model like all-MiniLM-L6-v2)

For activation profiles:
- **Per-layer divergence magnitude** (L2 norm of activation difference)
- **Per-layer relative divergence** (L2 of difference / L2 of activation)
- **Amplification factor** (ratio of layer k divergence to layer k-1 divergence)

### Expected Outcomes
- Control prompts: near-zero divergence (exact match ~99%+)
- Experimental prompts: measurable divergence, concentrated at low-margin positions
- Activation divergence should grow monotonically through layers, possibly with some layers amplifying more than others (especially attention layers)
- BF16/FP16 should show higher divergence rates than FP32

---

## Experiment 2: Training Divergence (Scenario 2)

### Goal
Measure how training the same model on different hardware produces weight-level and output-level divergence.

### Setup

**Fine-tuning configuration (identical across all platforms):**
- Base model: Pythia-410M (same checkpoint, verified identical)
- Method: LoRA (rank=8, alpha=16, applied to attention projection layers: q_proj, v_proj)
- Dataset: First 1,000 examples from the Alpaca dataset (tatsu-lab/alpaca), sorted deterministically by index
- Batch size: 4
- Learning rate: 2e-4
- Epochs: 3
- Optimizer: AdamW
- Random seed: 42 (set for Python, NumPy, PyTorch, and CUDA)
- Dataloader: shuffle=False (deterministic ordering)
- PyTorch deterministic mode: ON where supported (`torch.use_deterministic_algorithms(True)`)

The point of the aggressive determinism settings is to ensure that any divergence is attributable to hardware, not software randomness.

### Procedure

**Phase 2A — Training:**

On each platform (H_nvidia, H_mac, H_cloud):
1. Fine-tune the model using the exact configuration above
2. Save the resulting LoRA adapter weights
3. Record training loss at every step
4. Hash the final adapter weights (SHA-256)

Verify that the three adapter weight files have DIFFERENT hashes (they should, confirming hardware-induced training divergence).

**Phase 2B — Weight-level analysis:**

Compare the three sets of LoRA adapter weights pairwise:
- Element-wise absolute difference
- L2 distance between full adapter weight vectors
- Max absolute deviation
- Relative deviation (difference / weight magnitude)
- Per-layer weight divergence (do some layers diverge more than others?)

Compare training loss curves — they should be nearly identical but not bitwise equal.

**Phase 2C — Output-level divergence (the 2×2 matrix):**

Run all 400 prompts on all combinations:

| | Infer on H_nvidia | Infer on H_mac | Infer on H_cloud |
|---|---|---|---|
| **M_nvidia** (trained on nvidia) | Cell AA (baseline) | Cell AB | Cell AC |
| **M_mac** (trained on mac) | Cell BA | Cell BB (baseline) | Cell BC |
| **M_cloud** (trained on cloud) | Cell CA | Cell CB | Cell CC (baseline) |

This is a 3×3 matrix, but the key comparisons are:
- **Diagonal (AA, BB, CC):** Same training, same inference hardware. Run each twice for within-condition baseline.
- **Same row, different column:** Same training, different inference hardware — isolates inference divergence for a specific model.
- **Same column, different row:** Different training, same inference hardware — isolates training divergence.
- **Off-diagonal, different row AND column:** Both sources combined.

Record the same measurements as Experiment 1 for every cell.

**Phase 2D — Activation divergence for training-diverged models:**

For the 50-prompt activation subset, compare activation profiles between:
- M_nvidia vs M_mac, both running on H_nvidia (isolates training divergence at the activation level)
- M_nvidia on H_nvidia vs M_nvidia on H_mac (isolates inference divergence — should match Experiment 1)

This lets you see whether training-induced and inference-induced divergence affect different layers or have different amplification patterns.

### Expected Outcomes
- Weight-level divergence on the order of 10^-6 relative (based on Shanmugavelu et al.)
- Training divergence should produce LARGER output divergence than inference divergence alone
- The two sources should be approximately additive (but this is an empirical question)
- Loss curves should be visually indistinguishable despite weight-level differences

---

## Experiment 3: Divergence Prediction (Scenario 1, Follow-up)

### Goal
Determine whether layer-level activation divergence can predict output-level token divergence, and identify which layers are most informative.

### Prerequisites
Completed Experiments 1 and 2, specifically the activation recording phases.

### Procedure

1. From the activation data, construct a dataset:
   - For each prompt × each token position × each hardware pair:
     - Input features: per-layer activation divergence at that position (24-dimensional vector for Pythia-410M, one value per layer)
     - Label: did the output token diverge at this position? (binary)

2. Train a simple logistic regression on this dataset. This tells you: given the pattern of activation divergence across layers at a given position, can you predict whether the token will diverge?

3. Examine the learned coefficients. Large positive weights indicate layers where divergence is predictive of output divergence. This identifies your candidate intervention layers for fx(EC).

4. Compute precision/recall. High precision means: when the probe says "this will diverge," it usually does. High recall means: when divergence happens, the probe usually catches it.

### Additional Analysis

- Compute the correlation between per-layer divergence magnitude and output divergence probability. This is the simpler, non-learned version.
- Test whether a threshold rule works: "if divergence at layer k exceeds X, the output will diverge." Sweep thresholds to find the best layer and value.
- Check whether the predictive layer differs between attention and feedforward sublayers.

---

## Infrastructure and Code Notes

### Framework
- Use PyTorch throughout
- Use Hugging Face Transformers for model loading and generation
- Use PEFT library for LoRA fine-tuning
- Use `transformers.GenerationConfig` with `do_sample=False, temperature=0, top_k=1` for greedy decoding

### Activation Recording
- Use PyTorch forward hooks: `model.layers[i].register_forward_hook(hook_fn)` to capture each layer's output
- Store activations as FP32 tensors regardless of model precision
- Save to disk as `.pt` files (torch.save) — they compress well
- Expect ~50-100MB per prompt at 256 tokens for Pythia-410M with all layers recorded. For the 50-prompt subset, budget ~5GB of disk per platform.

### Data Storage and Organization
```
experiment_root/
  hardware_manifest.json
  prompts/
    candidates_2000.jsonl
    experimental_200.jsonl
    control_200.jsonl
    fragility_scores.csv
  experiment_0/
    logit_traces/
      prompt_0000.pt
      ...
  experiment_1/
    outputs/
      h_nvidia_run1.jsonl
      h_nvidia_run2.jsonl
      h_mac_run1.jsonl
      ...
    activations/
      h_nvidia/
        prompt_0000.pt
        ...
      h_mac/
        ...
    analysis/
      token_agreement.csv
      divergence_profiles.csv
      first_divergence_points.csv
  experiment_2/
    training/
      h_nvidia/
        adapter_weights/
        training_log.csv
      h_mac/
        ...
    weight_analysis/
      pairwise_distances.csv
      per_layer_divergence.csv
    outputs/
      (same structure as experiment_1 but for each M×H combination)
    analysis/
      ...
  experiment_3/
    probe_dataset.csv
    probe_results.json
```

### Reproducibility Checklist
Before starting any experiment run:
- [ ] Hardware manifest recorded
- [ ] Model checkpoint hashes verified (must match across platforms)
- [ ] PyTorch version identical across platforms (or documented if not possible)
- [ ] Random seeds set: `torch.manual_seed(42)`, `np.random.seed(42)`, `random.seed(42)`
- [ ] CUDA deterministic mode: `torch.use_deterministic_algorithms(True)` (on NVIDIA platforms)
- [ ] CUBLAS workspace config: `os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'` (on NVIDIA)
- [ ] Dataloader workers: 0 (to avoid multiprocessing nondeterminism)
- [ ] Prompts loaded from the fixed .jsonl files (not regenerated)

### Rough Time Budget

| Task | H_nvidia | H_mac | H_cloud |
|------|----------|-------|---------|
| Experiment 0 (prompt scoring, 2000 prompts) | ~3 hrs | — | — |
| Experiment 1A (6 inference runs × 400 prompts × 256 tokens) | ~4 hrs | ~5 hrs | ~2 hrs |
| Experiment 1B (activation recording, 50 prompts) | ~2 hrs | ~3 hrs | ~1 hr |
| Experiment 1C (repeat at BF16/FP16) | ~4 hrs | ~5 hrs | ~2 hrs |
| Experiment 2A (fine-tuning, 3 epochs on 1000 examples) | ~30 min | ~45 min | ~15 min |
| Experiment 2C (9+ inference runs) | ~6 hrs | ~8 hrs | ~3 hrs |
| Experiment 2D (activation recording) | ~2 hrs | ~3 hrs | ~1 hr |
| Experiment 3 (analysis, probe training) | ~30 min | — | — |

Total wall-clock with parallelism across machines: approximately 1.5-2 days.
Cloud cost estimate: 5-8 hours of A100 time ≈ $10-20.

### If You Need to Reduce Scope
Priority order (do these first if time is tight):
1. Experiment 0 (prompt selection) — necessary for everything else
2. Experiment 1A at FP32 (core inference divergence measurement)
3. Experiment 1B (activation profiles — the most novel data)
4. Experiment 2A + 2B (training weight divergence)
5. Experiment 2C (output-level training divergence)
6. Everything else

You can skip the cloud platform entirely if budget is a concern — the nvidia-vs-mac comparison is the most architecturally interesting pair anyway.

---

## Analysis Deliverables

After all experiments, the key outputs are:

1. **Divergence rate table:** For each hardware pair, at each precision, what fraction of prompts produce identical outputs? Broken down by prompt category.
2. **Divergence onset distribution:** Histogram of first-divergence token positions. Does divergence tend to happen early or late in generation?
3. **Margin-divergence correlation:** Plot of logit margin at divergence point vs whether divergence occurred. Expected to show a clear threshold.
4. **Layer divergence profiles:** Curves showing activation divergence magnitude vs layer index, averaged across prompts, for each hardware pair. This is the most visually informative result.
5. **Amplification factors:** Which layers amplify divergence most? Attention vs feedforward?
6. **Training vs inference variance decomposition:** From the 2×2 design, how much of total output variance comes from training vs inference hardware differences?
7. **Probe results:** Can layer-level divergence predict output divergence? If so, which layers are most predictive?
