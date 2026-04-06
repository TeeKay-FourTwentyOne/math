"""
Experiment C: Micro-benchmark of FP reduction order.

Strip away the model. Take known vectors, sum them with different
chunk sizes, and measure how the result changes. This is the atomic
unit of all cross-platform divergence.

Also test: matrix multiplication with different algorithms.
"""

import torch
import numpy as np
import os
import random

SEED = 42

def setup_determinism():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)

setup_determinism()
device = torch.device("cuda")

print("=" * 70)
print("MICRO-BENCHMARK: FP Reduction Order Effects")
print("=" * 70)

# ============================================================
# TEST 1: Summation with different chunk sizes
# ============================================================
print("\n--- TEST 1: Vector summation with different chunk sizes ---")
print("Sum 4096 random FP32 values using different groupings.\n")

# Create a known vector (same on all platforms via seed)
torch.manual_seed(42)
vec = torch.randn(4096, dtype=torch.float32, device="cpu")  # CPU for control

# "True" sum (sequential, CPU, FP64 for reference)
true_sum = vec.double().sum().item()
print(f"Reference sum (FP64 sequential): {true_sum:.15f}")

# Sequential sum (CPU FP32)
seq_sum = 0.0
for v in vec.tolist():
    seq_sum += v
print(f"Sequential FP32 (Python loop):   {seq_sum:.15f}  diff={abs(seq_sum - true_sum):.2e}")

# PyTorch CPU sum
cpu_sum = vec.sum().item()
print(f"PyTorch CPU sum:                 {cpu_sum:.15f}  diff={abs(cpu_sum - true_sum):.2e}")

# PyTorch CUDA sum (deterministic mode)
cuda_vec = vec.to(device)
cuda_sum_det = cuda_vec.sum().item()
print(f"PyTorch CUDA sum (deterministic): {cuda_sum_det:.15f}  diff={abs(cuda_sum_det - true_sum):.2e}")

# Now chunk manually and compare
print(f"\nChunked summation (sum chunks of size K, then sum the partial sums):")
print(f"{'Chunk K':>8s} {'Result':>20s} {'Diff from ref':>15s} {'Diff from CUDA':>15s}")
for chunk_size in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]:
    # Chunk the vector, sum each chunk, then sum the partial sums
    chunks = vec.reshape(-1, chunk_size) if 4096 % chunk_size == 0 else vec[:4096 // chunk_size * chunk_size].reshape(-1, chunk_size)
    partial_sums = chunks.sum(dim=1)
    chunked_result = partial_sums.sum().item()
    diff_ref = abs(chunked_result - true_sum)
    diff_cuda = abs(chunked_result - cuda_sum_det)
    print(f"{chunk_size:8d} {chunked_result:20.15f} {diff_ref:15.2e} {diff_cuda:15.2e}")

# ============================================================
# TEST 2: Same vector, CUDA vs CPU reduction
# ============================================================
print("\n--- TEST 2: CUDA vs CPU reduction on same data ---")
print("Different vector sizes to see how error scales.\n")

print(f"{'Size':>8s} {'CPU sum':>20s} {'CUDA sum':>20s} {'|Diff|':>12s} {'Rel diff':>12s}")
for size in [64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 65536]:
    torch.manual_seed(42)
    v = torch.randn(size, dtype=torch.float32)
    cpu_s = v.sum().item()
    cuda_s = v.to(device).sum().item()
    diff = abs(cpu_s - cuda_s)
    rel = diff / abs(cpu_s) if cpu_s != 0 else 0
    print(f"{size:8d} {cpu_s:20.10f} {cuda_s:20.10f} {diff:12.2e} {rel:12.2e}")

# ============================================================
# TEST 3: Matrix multiplication — CPU vs CUDA
# ============================================================
print("\n--- TEST 3: Matrix multiplication CPU vs CUDA ---")
print("Simulates a single transformer layer's computation.\n")

print(f"{'M':>5s} {'K':>5s} {'N':>5s} {'Max |diff|':>12s} {'Mean |diff|':>12s} {'Rel max':>12s}")

# Test dimensions matching Pythia-410M
test_cases = [
    (1, 1024, 1024, "hidden @ hidden (attention)"),
    (1, 1024, 4096, "hidden @ up_proj (MLP)"),
    (1, 4096, 1024, "mlp_out @ down_proj"),
    (1, 1024, 50304, "hidden @ lm_head (logits!)"),
    (256, 1024, 1024, "full seq @ hidden"),
    (256, 1024, 50304, "full seq @ lm_head"),
]

for M, K, N, desc in test_cases:
    torch.manual_seed(42)
    A = torch.randn(M, K, dtype=torch.float32)
    B = torch.randn(K, N, dtype=torch.float32)

    cpu_result = A @ B
    cuda_result = (A.to(device) @ B.to(device)).cpu()

    diff = (cpu_result - cuda_result).abs()
    max_diff = diff.max().item()
    mean_diff = diff.mean().item()
    rel_max = max_diff / cpu_result.abs().max().item() if cpu_result.abs().max().item() > 0 else 0

    print(f"{M:5d} {K:5d} {N:5d} {max_diff:12.2e} {mean_diff:12.2e} {rel_max:12.2e}  {desc}")

# ============================================================
# TEST 4: The actual logit computation path
# ============================================================
print("\n--- TEST 4: Simulated logit computation (layernorm -> matmul -> top-k) ---")
print("Does the top-1 token ever change between CPU and CUDA?\n")

torch.manual_seed(42)
# Simulate: hidden state -> layer_norm -> lm_head -> logits -> argmax
hidden = torch.randn(1, 256, 1024, dtype=torch.float32)
ln_weight = torch.randn(1024, dtype=torch.float32)
ln_bias = torch.randn(1024, dtype=torch.float32)
lm_head = torch.randn(50304, 1024, dtype=torch.float32)

# CPU path
h_cpu = torch.nn.functional.layer_norm(hidden, [1024], ln_weight, ln_bias)
logits_cpu = h_cpu @ lm_head.t()  # (1, 256, 50304)
top1_cpu = logits_cpu.argmax(dim=-1)  # (1, 256)

# CUDA path
h_cuda = torch.nn.functional.layer_norm(
    hidden.to(device), [1024], ln_weight.to(device), ln_bias.to(device))
logits_cuda = h_cuda @ lm_head.to(device).t()
top1_cuda = logits_cuda.argmax(dim=-1).cpu()

# Compare
logit_diff = (logits_cpu - logits_cuda.cpu()).abs()
token_mismatches = (top1_cpu != top1_cuda).sum().item()

print(f"Logit differences across 256 positions x 50304 vocab:")
print(f"  Max:  {logit_diff.max().item():.6f}")
print(f"  Mean: {logit_diff.mean().item():.6f}")
print(f"  Top-1 token mismatches: {token_mismatches}/256")

# For mismatched positions, show the margin
if token_mismatches > 0:
    mismatch_positions = (top1_cpu != top1_cuda).squeeze().nonzero().squeeze().tolist()
    if isinstance(mismatch_positions, int):
        mismatch_positions = [mismatch_positions]
    print(f"\n  Mismatched positions:")
    for pos in mismatch_positions[:10]:
        cpu_logits = logits_cpu[0, pos]
        cuda_logits = logits_cuda[0, pos].cpu()
        cpu_top2 = cpu_logits.topk(2)
        cuda_top2 = cuda_logits.topk(2)
        margin_cpu = (cpu_top2.values[0] - cpu_top2.values[1]).item()
        margin_cuda = (cuda_top2.values[0] - cuda_top2.values[1]).item()
        print(f"    pos {pos}: CPU tok={top1_cpu[0,pos].item()} CUDA tok={top1_cuda[0,pos].item()} "
              f"CPU margin={margin_cpu:.6f} CUDA margin={margin_cuda:.6f}")
else:
    # Show how close the closest call was
    margins_cpu = logits_cpu.topk(2, dim=-1).values
    margins_diff = (margins_cpu[0,:,0] - margins_cpu[0,:,1])
    min_margin_pos = margins_diff.argmin().item()
    min_margin = margins_diff[min_margin_pos].item()
    logit_diff_at_min = logit_diff[0, min_margin_pos].max().item()
    print(f"\n  All tokens match. Closest call at position {min_margin_pos}:")
    print(f"    Margin: {min_margin:.6f}, Max logit diff at that position: {logit_diff_at_min:.6f}")
    print(f"    Safety factor: {min_margin / logit_diff_at_min:.1f}x")
