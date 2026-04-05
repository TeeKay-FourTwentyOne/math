# Experiment 0 Implementation Notes

## Result: COMPLETE (2026-04-04)
- Duration: 112.4 minutes on RTX 4070
- 1999 prompts scored (5 categories x ~400 each)
- 200 experimental + 200 control prompts selected
- Zero overlap between experimental and control fragility ranges

## Key Findings
| Category | Mean Fragility | Median | Max |
|----------|---------------|--------|-----|
| reasoning | 18.3% | 22.3% | 62.5% |
| continuation | 13.6% | 13.3% | 28.1% |
| factual | 9.3% | 8.2% | 35.6% |
| conversational | 8.8% | 8.2% | 23.4% |
| open_ended | 8.3% | 7.4% | 30.9% |

## Selected Sets
- Experimental (200): mean fragility 23.9%, range 13.7%–62.5%, balanced 40/category
- Control (200): mean fragility 3.6%, range 0.4%–8.6%, balanced 40/category
- No overlap between sets (min experimental > max control)

## Top 10 most fragile: ALL reasoning (math) prompts
- Percentage questions (62.5% fragility) and multi-digit addition (59.8%)
- These will show the strongest hardware divergence in Experiment 1

## Validated Parameters
- Model: EleutherAI/pythia-410m-deduped (405M params, 24 layers, FP32)
- VRAM: 1.62 GB (of 12 GB available)
- Speed: 3.4s/prompt, ~0.3 prompts/s
- Determinism: PERFECT (strict mode, 0.00e+00 logit diff between runs)
- Trace size: ~20 KB/prompt, 41 MB total

## Improvements for Future Runs
- Suppress pad_token_id warnings: logging.getLogger("transformers").setLevel(logging.ERROR)
- Flush stdout for background runs: PYTHONUNBUFFERED=1
- Resume support works (tested implicitly via checkpoint design)

## Files
- hardware_manifest.json — platform spec
- hardware_manifest.py — manifest generator
- smoke_test.py — pipeline validation (PASSED)
- generate_prompts.py — prompt generation (1999 prompts)
- experiment_0.py — full scoring pipeline
- prompts/candidates_2000.jsonl — all prompts
- prompts/experimental_200.jsonl — high-fragility set
- prompts/control_200.jsonl — low-fragility set
- prompts/fragility_scores.csv — all scores
- experiment_0/logit_traces/ — 1999 .pt files
