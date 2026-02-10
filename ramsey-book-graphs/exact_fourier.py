"""
Exact Fourier characterization: for the "standard" regime (n >= 12),
verify the exact values of P(0)/m and the constraint margins.

Key finding from fourier_constraints.py:
- For n >= 12: avg A(d)+B(d) = n exactly (integer!)
- Red threshold = n-2, so excess = 2 above threshold on average
- Red fluctuation needed: exactly -P(0)/m + (n-2) = -(n + 1/(2n-1) - 1 + ...) + n - 2 ~ -1.5
- The constraint is VERY tight: red_slack / |D11| is small

Let's compute the exact average and verify it equals n.
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta

KNOWN = {
    12: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18}, "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}},
    14: {"D11": {5, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20, 22}, "D12": {0, 1, 2, 7, 8, 10, 13, 14, 17, 18, 21, 23, 25}},
    16: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25}, "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}},
    18: {"D11": {6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20, 21, 22, 24, 25, 26, 27, 29}, "D12": {0, 1, 2, 3, 8, 9, 10, 12, 14, 15, 19, 20, 22, 24, 25, 28, 32}},
    20: {"D11": {7, 8, 9, 10, 12, 13, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27, 29, 30, 31, 32}, "D12": {0, 1, 2, 5, 6, 8, 12, 14, 15, 17, 18, 22, 23, 25, 27, 28, 33, 36, 37}},
}

n22_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solution_n22.json")
if os.path.exists(n22_path):
    with open(n22_path) as f:
        data = json.load(f)
    KNOWN[22] = {"D11": set(data["parameters"]["D11"]), "D12": set(data["parameters"]["D12"])}


print("EXACT AVERAGE VERIFICATION")
print("=" * 80)

for n in sorted(KNOWN.keys()):
    m = 2 * n - 1
    D11 = KNOWN[n]["D11"]
    D12 = KNOWN[n]["D12"]

    # sum_{d=1}^{m-1} (A(d)+B(d)) = |D11|(|D11|-1) + |D12|(|D12|-1)
    total = len(D11) * (len(D11) - 1) + len(D12) * (len(D12) - 1)
    avg = total / (m - 1)

    # For |D11| = n, |D12| = n-1:
    # total = n(n-1) + (n-1)(n-2) = (n-1)(n + n - 2) = (n-1)(2n-2) = 2(n-1)^2
    # avg = 2(n-1)^2 / (m-1) = 2(n-1)^2 / (2n-2) = 2(n-1)^2 / (2(n-1)) = n-1
    # Wait, that gives n-1, not n!
    theoretical = 2 * (n - 1) ** 2 / (2 * (n - 1))

    # Hmm but the output showed avg = n for n >= 12. Let me recheck.
    # Oh wait: the output from fourier_constraints.py said "Average A(d)+B(d) over d=1..m-1: 11.000" for n=12.
    # That's computed as total/(m-1).
    # For n=12: |D11|=12, |D12|=11
    # total = 12*11 + 11*10 = 132 + 110 = 242
    # avg = 242/22 = 11.0 = n-1. Not n! It's n-1!

    print(f"n={n}: |D11|={len(D11)}, |D12|={len(D12)}, m={m}")
    print(f"  total = {total}, avg = {avg:.4f}, n-1 = {n-1}, n = {n}")
    print(f"  theoretical (2(n-1)^2/(2(n-1))) = {theoretical:.4f}")

    # Wait, I was confused. Let me re-examine. The average over d=1..m-1:
    # But Delta(S,S,0) = |S|, and we're summing d=1..m-1.
    # sum_{d=0}^{m-1} Delta(S,S,d) = |S|^2 (all ordered pairs)
    # sum_{d=1}^{m-1} Delta(S,S,d) = |S|^2 - |S| = |S|(|S|-1)
    # So avg_{d=1..m-1} Delta(S,S,d) = |S|(|S|-1)/(m-1)
    #
    # For D11: |D11|(|D11|-1)/(m-1) = n(n-1)/(2n-2) = n/2
    # For D12: |D12|(|D12|-1)/(m-1) = (n-1)(n-2)/(2n-2) = (n-2)/2
    # Total: n/2 + (n-2)/2 = (2n-2)/2 = n-1.
    # So avg(A+B) = n-1 over {1,...,m-1}. The output says n-1 too (11 for n=12).

    # But the average of A+B restricted to D11 (red edges):
    red_sum = sum(Delta(D11, D11, d, m) + Delta(D12, D12, d, m) for d in D11)
    red_avg = red_sum / len(D11)

    # And restricted to complement (blue edges):
    blue_diffs = [d for d in range(1, m) if d not in D11]
    blue_sum = sum(Delta(D11, D11, d, m) + Delta(D12, D12, d, m) for d in blue_diffs)
    blue_avg = blue_sum / len(blue_diffs)

    print(f"  Red avg (d in D11): {red_avg:.4f}, threshold: {n-2}")
    print(f"  Blue avg (d not in D11): {blue_avg:.4f}, blue max threshold: {n-1 + 2*(len(D11)+len(D12)) - (4*n-4)}")
    print(f"  Overall avg: {avg:.4f}")
    print(f"  Red avg - overall avg = {red_avg - avg:.4f}")
    print(f"  Red avg - threshold = {red_avg - (n-2):.4f}")

    # Key: red_avg - threshold is the average "excess" of red common neighbors
    # For n=12: red_avg = 118/12 = 9.833, threshold = 10. So excess = -0.167 (below threshold!)
    # Oh wait, threshold for red is n-2 = 10, and red_avg < 10. So on average, red is fine.
    # The constraint is about the MAX over d in D11, not the average.

    # Distribution of A(d)+B(d) for d in D11:
    red_vals = sorted([Delta(D11, D11, d, m) + Delta(D12, D12, d, m) for d in D11])
    print(f"  Red A+B values: {red_vals}")
    print(f"  Max red = {max(red_vals)}, threshold = {n-2}")

print("\n" + "=" * 80)
print("EXACT CONSTRAINT TIGHTNESS")
print("=" * 80)
print("""
For the standard regime (n >= 12, |D11| = n, |D12| = n-1):

P(0)/m = (n^2 + (n-1)^2) / (2n-1) = (2n^2 - 2n + 1)/(2n-1)
       = n - 1/(2(2n-1)) ... let me compute exactly.

(2n^2 - 2n + 1)/(2n-1): do polynomial division.
2n^2 - 2n + 1 = (2n-1) * n + (-2n + 1 + n) = (2n-1)*n + (-n + 1)
Hmm: (2n-1)*n = 2n^2 - n. Remainder: -2n+1 - (-n) = -n + 1.
So P(0)/m = n + (-n+1)/(2n-1) = n - (n-1)/(2n-1).

For large n: P(0)/m ~ n - 1/2.

Red fluctuation bound: F(d) <= (n-2) - P(0)/m = (n-2) - n + (n-1)/(2n-1) = -2 + (n-1)/(2n-1)
For large n: ~ -2 + 1/2 = -3/2 = -1.5.

So we need: for all d in D11, F(d) <= -3/2 + O(1/n).

Blue fluctuation bound: F(d) <= (n+1) - P(0)/m = 1 + (n-1)/(2n-1)
For large n: ~ 1 + 1/2 = 3/2.

So we need: for all d not in D11, F(d) <= 3/2 + O(1/n).

The fluctuation F(d) has mean 0 (over all d != 0) and the red/blue bounds are
approximately ±3/2. The std of F is ~1.5-2 based on the data. So the constraint
asks that ~half the values are below -3/2 and the other half are below +3/2.
This is a 0.5-1 sigma event on each side -- quite feasible!

The real question for a proof: can we show that for EVERY symmetric D11 with
|D11| = (m+1)/2 (where m = prime ≡ 3 mod 4), there exists a suitable D12?
Or equivalently: given the power spectrum |hat{D11}(k)|^2, can we always find
D12 with |D12| = (m-1)/2 whose power spectrum |hat{D12}(k)|^2 "corrects" the
fluctuations to satisfy the constraint?
""")
