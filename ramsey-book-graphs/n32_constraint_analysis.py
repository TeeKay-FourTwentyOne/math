"""
Deep constraint analysis for n=32, m=63.

Investigate why |D11|=32 and |D11|=30 don't work.
The V1V1 constraint for red edge at d: Delta(D11,D11,d) + Delta(D12,D12,d) <= 30
The V1V1 constraint for blue edge at d: (N-2) - 2*d1 + Delta(D11,D11,d) + Delta(D12,D12,d) <= 31

Since red_thresh = 30 and blue_thresh = 31, and N-2 = 124, d1 = |D11| + |D12|:
- Red: A(d) + B(d) <= 30  where A=Delta(D11,D11,d), B=Delta(D12,D12,d)
- Blue: 124 - 2*d1 + A(d) + B(d) <= 31  =>  A(d) + B(d) <= 2*d1 - 93

For |D11|=32, |D12|=31: d1 = 63, so blue constraint: A(d)+B(d) <= 33
For |D11|=30, |D12|=31: d1 = 61, so blue constraint: A(d)+B(d) <= 29

Wait - for |D11|=30, the blue constraint on V1V1 is TIGHTER than the red constraint!
Red: A(d)+B(d) <= 30
Blue: A(d)+B(d) <= 29

This means for |D11|=30, blue V1V1 edges need even lower autocorrelation!

Let me also work out V2V2:
V2V2(d) = A(d) + B(-d) + (m-2-2|D11|) + 2[d in D11]

For |D11|=32: V2V2_const = 63-2-64 = -3
  V2V2 red (d not in D11): A(d) + B(-d) - 3 <= 30  =>  A(d)+B(-d) <= 33
  V2V2 blue (d in D11): blue = (N-2)-2*d2 + V2V2(d)
    d2 = |D22|+|D12| = 30+31 = 61
    blue = 124-122 + A(d)+B(-d)-3+2 = 1 + A(d)+B(-d)
    Constraint: 1 + A(d)+B(-d) <= 31  =>  A(d)+B(-d) <= 30

For |D11|=30: V2V2_const = 63-2-60 = 1
  V2V2 red (d not in D11): A(d) + B(-d) + 1 <= 30  =>  A(d)+B(-d) <= 29
  V2V2 blue (d in D11): blue = (N-2)-2*d2 + V2V2(d)
    d2 = |D22|+|D12| = 32+31 = 63
    blue = 124-126 + A(d)+B(-d)+1+2 = 1 + A(d)+B(-d)
    Constraint: 1 + A(d)+B(-d) <= 31  =>  A(d)+B(-d) <= 30

Summary of constraints for |D11|=32, |D12|=31:
  V1V1 red (d in D11):   A(d)+B(d) <= 30
  V1V1 blue (d not in D11): A(d)+B(d) <= 33
  V2V2 blue (d in D11):  A(d)+B(m-d) <= 30
  V2V2 red (d not in D11): A(d)+B(m-d) <= 33

So the binding constraints are:
  For d in D11: A(d)+B(d) <= 30 AND A(d)+B(m-d) <= 30
  For d not in D11: (non-binding, just <= 33)

Summary for |D11|=30, |D12|=31:
  V1V1 red (d in D11):   A(d)+B(d) <= 30
  V1V1 blue (d not in D11): A(d)+B(d) <= 29  [TIGHTER!]
  V2V2 blue (d in D11):  A(d)+B(m-d) <= 30
  V2V2 red (d not in D11): A(d)+B(m-d) <= 29  [TIGHTER!]

For |D11|=30, the BLUE constraints are the binding ones!

Let's compute the sum constraint. Sum of A(d) over all d = |D11|^2.
Sum of B(d) over all d = |D12|^2.

For |D11|=32: sum A = 1024, average A = 1024/63 ≈ 16.25 (for d=1..62)
  A(0) = 32, so sum(A[1:]) = 992, average A(d) for d!=0 = 992/62 ≈ 16.0
For |D12|=31: sum B = 961, B(0) = 31, sum(B[1:]) = 930, average B(d) for d!=0 = 930/62 ≈ 15.0
Average A(d)+B(d) for d!=0 ≈ 31.0

The threshold for d in D11 is 30. So on average A(d)+B(d) ≈ 31 but we need it <= 30!
This means the average exceeds the threshold by ~1 for red edges.

This is the fundamental issue: the average autocorrelation sum is just barely above the threshold.

For |D11|=30:
  sum A = 900, A(0)=30, sum(A[1:]) = 870, avg A(d) for d!=0 = 870/62 ≈ 14.03
  sum B = 961, B(0)=31, sum(B[1:]) = 930, avg B(d) for d!=0 = 930/62 ≈ 15.0
  Average A(d)+B(d) ≈ 29.03
  Blue threshold (d not in D11): 29
  Average is right at the threshold for blue edges too!
"""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

m = 63
n = 32

print("=== Constraint Analysis for n=32, m=63 ===\n")

for d11_size in [32, 30, 28, 34]:
    d12_size = 31
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m

    print(f"\n--- |D11|={d11_size}, |D12|={d12_size}, |D22|={d22_size} ---")
    print(f"d1={d1}, d2={d2}, N={N}")

    # V1V1 constraints
    v11_red_thresh_sum = n - 2  # 30
    v11_blue_rhs = 2 * d1 - (N - 2) + (n - 1)  # A+B <= this
    print(f"V1V1 red (d in D11):    A(d)+B(d) <= {v11_red_thresh_sum}")
    print(f"V1V1 blue (d not D11):  A(d)+B(d) <= {v11_blue_rhs}")

    # V2V2 constraints
    v22_const = m - 2 - 2 * d11_size
    # V2V2 red: d not in D11 (d in D22)
    # V2V2(d) = A(d) + B(m-d) + v22_const
    # A(d) + B(m-d) + v22_const <= 30
    v22_red_rhs = n - 2 - v22_const
    # V2V2 blue: d in D11
    # blue = (N-2) - 2*d2 + A(d) + B(m-d) + v22_const + 2
    # = (N-2) - 2*d2 + v22_const + 2 + A(d) + B(m-d) <= n-1
    v22_blue_rhs = (n - 1) - ((N - 2) - 2 * d2 + v22_const + 2)
    print(f"V2V2 red (d not D11):   A(d)+B(m-d) <= {v22_red_rhs}")
    print(f"V2V2 blue (d in D11):   A(d)+B(m-d) <= {v22_blue_rhs}")

    # Average autocorrelation sums
    avg_A = (d11_size ** 2 - d11_size) / (m - 1)
    avg_B = (d12_size ** 2 - d12_size) / (m - 1)
    avg_sum = avg_A + avg_B

    print(f"\nAverage A(d) for d!=0: {avg_A:.2f}")
    print(f"Average B(d) for d!=0: {avg_B:.2f}")
    print(f"Average A(d)+B(d): {avg_sum:.2f}")

    # Binding constraint analysis
    binding = min(v11_red_thresh_sum, v11_blue_rhs, v22_red_rhs, v22_blue_rhs)
    print(f"\nTightest constraint on A(d)+B(?): {binding}")
    print(f"Average vs tightest: {avg_sum:.2f} vs {binding}")
    slack = binding - avg_sum
    print(f"Average slack: {slack:.2f} ({'feasible on average' if slack >= 0 else 'INFEASIBLE on average'})")

    # For d in D11: need A(d)+B(d) <= min(v11_red, ...) and A(d)+B(m-d) <= min(v22_blue, ...)
    # For d not in D11: need A(d)+B(d) <= min(v11_blue, ...) and A(d)+B(m-d) <= min(v22_red, ...)
    d11_binding = min(v11_red_thresh_sum, v22_blue_rhs)
    d22_binding = min(v11_blue_rhs, v22_red_rhs)
    print(f"\nFor d in D11 (|D11|={d11_size} differences):")
    print(f"  Need A(d)+B(d) <= {v11_red_thresh_sum} AND A(d)+B(m-d) <= {v22_blue_rhs}")
    print(f"  Expected: {d11_size} constraints at threshold {d11_binding}")
    print(f"For d not in D11 (|D22|={d22_size} differences):")
    print(f"  Need A(d)+B(d) <= {v11_blue_rhs} AND A(d)+B(m-d) <= {v22_red_rhs}")
    print(f"  Expected: {d22_size} constraints at threshold {d22_binding}")

    # Total budget: sum of all A(d)+B(d) for d=1..m-1 is fixed
    total_sum = (d11_size ** 2 - d11_size) + (d12_size ** 2 - d12_size)
    print(f"\nTotal sum(A(d)+B(d), d=1..{m-1}) = {total_sum}")
    # If all d in D11 achieve exactly the threshold:
    best_case_d11 = d11_size * d11_binding
    remaining = total_sum - best_case_d11
    per_d22 = remaining / d22_size if d22_size > 0 else 0
    print(f"If all d in D11 at threshold {d11_binding}: budget used = {best_case_d11}")
    print(f"Remaining for d not in D11: {remaining}")
    print(f"Average per d22: {per_d22:.2f} (threshold = {d22_binding})")

    # Also check the cross-correlation budget
    # A(d)+B(m-d) has a different average since B(m-d) has same distribution as B(d)
    # Total sum of A(d)+B(m-d) = sum(A(d)) + sum(B(m-d)) = same as sum(A(d)) + sum(B(d))
    print(f"\nTotal sum(A(d)+B(m-d), d=1..{m-1}) = {total_sum} (same)")

print("\n\n=== Key Insight ===")
print("For |D11|=32: average A(d)+B(d) ≈ 31, but red threshold = 30.")
print("Need ALL 32 red-edge differences to be at most 30, when average is 31.")
print("This requires high variance / lucky distribution -- hard but not impossible.")
print()
print("For |D11|=30: average A(d)+B(d) ≈ 29, but blue threshold = 29.")
print("Need ALL 32 blue-edge differences to be at most 29, when average is 29.")
print("Same issue from the blue side.")
print()
print("Neither size has much slack. This explains why n=32 is hard.")
print("The 'sweet spot' where average is well below threshold doesn't exist for m=63.")

if __name__ == "__main__":
    pass
