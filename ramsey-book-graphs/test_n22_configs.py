"""Test initial configurations for n=22 search."""

from ramsey_core import BlockCirculantGraph, verify_construction, compute_violation_cost

def quadratic_residues(m):
    return {(x * x) % m for x in range(1, m)} - {0}

n = 22
m = 43
half_m = m // 2  # 21

QR43 = quadratic_residues(43)
NQR43 = set(range(1, 43)) - QR43

print(f"n={n}, m={m}, N={4*n-2}")
print(f"QR43 = {sorted(QR43)}")
print(f"|QR43| = {len(QR43)}")
print()

# Test different D11 seeds
D11_seeds = [
    ("D11={1..11}", set(range(1, 12))),
    ("D11=odd", set(range(1, 22, 2))),
    ("D11=QR43[:11]", {1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17}),
]

D12_options = [
    ("D12=QR43+{0}", QR43 | {0}),
    ("D12=QR43", QR43.copy()),
    ("D12=NQR43+{0}", NQR43 | {0}),
]

print("Testing configurations (D22 = complement of D11):")
print("=" * 60)

best_cost = float('inf')
best_config = None

for d11_name, D11_half in D11_seeds:
    for d12_name, D12 in D12_options:
        # D22_half = complement of D11_half in {1..21}
        D22_half = set(range(1, half_m + 1)) - D11_half

        G = BlockCirculantGraph(n=n, D11=D11_half, D12=D12, D22=D22_half)
        cost = compute_violation_cost(G)

        config_name = f"{d11_name}, {d12_name}"
        print(f"{config_name}:")
        print(f"  |D11|={len(G.D11)}, |D12|={len(G.D12)}, |D22|={len(G.D22)}")
        print(f"  Cost: {cost}")

        if cost < best_cost:
            best_cost = cost
            best_config = config_name

print()
print(f"Best starting config: {best_config} with cost {best_cost}")
