# Implementation Plan: Step 7 GCP Deployment for 6D Hénon-Heiles Benchmark

## Summary

Create a new `step7_gcp/` directory in `C:\Projects\physics` with self-contained scripts to deploy the LP+PINN solver to Google Cloud Platform for the 6D Hénon-Heiles benchmark.

**Constraints:**
- ADD only - no modifications to existing files (20+ hour training run in progress)
- Self-contained training script - no imports from existing project code
- All computation on GCP, not locally

---

## Files to Create

All files go in `C:\Projects\physics\step7_gcp\`:

| File | Purpose |
|------|---------|
| `README.md` | User instructions for GCP setup and deployment |
| `gcp_setup.sh` | Creates GCS bucket and Spot VM with L4 GPU |
| `startup_script.sh` | VM startup script - downloads code, starts training |
| `train_henon_heiles.py` | Self-contained training script (~600 lines) |
| `monitor.sh` | Check training progress from local machine |
| `cleanup.sh` | Tear down GCP resources when done |

---

## Key Implementation Details

### 1. `train_henon_heiles.py` Architecture

Based on the instructions template, with hyperparameters aligned to Step 6:

```
Network: 8 hidden layers x 512 neurons (matches Step 6)
Activation: Tanh
Input: 6 (position coordinates)
Output: 12 (real + imaginary momentum components)

Loss functions:
- physics_loss: QHJE residual (p² + iℏ∇·p = 2m(E-V))
- curl_loss: Irrotationality (∂pᵢ/∂xⱼ = ∂pⱼ/∂xᵢ)

Hyperparameters:
- Learning rate: 1e-4
- Collocation points: 100,000 per epoch
- Max epochs: 50,000
- Checkpoint interval: 1000 epochs
```

### 2. Hénon-Heiles Potential

```python
V(x₁,...,x₆) = ½Σᵢxᵢ² + λΣᵢ(xᵢ²xᵢ₊₁ - xᵢ₊₁³/3)
```

With periodic boundary (x₇ = x₁) and λ = 0.111803 (1/√80).

### 3. GCP Configuration

```
VM: g2-standard-8 (L4 GPU) or n1-standard-8 + T4
Provisioning: Spot (preemptible) - ~$0.25-0.40/hr
Region: us-central1-a
Boot disk: 100GB SSD with PyTorch Deep Learning VM image
```

### 4. Preemption Handling

- Checkpoints saved to GCS every 30 minutes via cron
- Startup script checks for existing checkpoints and resumes
- VM set to STOP (not DELETE) on termination

---

## Design Decision: Template vs Full Step 6 Architecture

The instructions provide a simplified template without:
- Learnable poles
- Supervision loss
- Quantization loss

**Recommendation: Use the template approach** because:
1. Hénon-Heiles has no obvious supervision target (unlike coupled oscillator with alpha_matrix)
2. Simpler architecture is easier to debug on GCP
3. Template is explicitly designed for this benchmark
4. Can always iterate if 2% accuracy target isn't met

---

## Verification

After deployment:

1. **VM health**: `./monitor.sh status` returns RUNNING
2. **Training started**: `./monitor.sh log` shows epoch progress
3. **Checkpoints**: `./monitor.sh checkpoints` shows files in GCS
4. **Energy convergence**: E should decrease from 3.0 toward ~2.97
5. **Final results**: Energy within 2% of 2.97 a.u. (target: 2.91-3.03)

**Success criteria** (from instructions):
- Final energy: 2.95-2.99 a.u.
- Below harmonic ZPE (3.0): YES
- Error vs ~2.97: < 2%

---

## Implementation Steps

1. Create `C:\Projects\physics\step7_gcp\` directory
2. Write `README.md` with prerequisites and workflow
3. Write `gcp_setup.sh` (bucket creation, VM creation)
4. Write `startup_script.sh` (download code, start training, cron sync)
5. Write `train_henon_heiles.py` (self-contained ~600 lines)
6. Write `monitor.sh` (log, checkpoints, results, ssh, status commands)
7. Write `cleanup.sh` (delete VM and bucket with confirmation)

---

## Estimated GCP Cost

| Configuration | Hourly (Spot) | 24hr Total |
|---------------|---------------|------------|
| g2-standard-8 (L4) | $0.25-0.40 | $6-10 |
| n1-standard-8 + T4 | $0.15-0.25 | $4-6 |

Expected training time: ~20 hours → **$5-8 total**
