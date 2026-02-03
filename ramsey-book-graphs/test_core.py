"""
Unit tests for ramsey_core.py

Verifies fast formula results against brute-force adjacency matrix computations.
"""

import unittest
import random
from ramsey_core import (
    BlockCirculantGraph, Delta, Sigma,
    common_neighbors_V1V1, common_neighbors_V2V2, common_neighbors_V1V2,
    common_neighbors_brute_force, verify_construction, make_symmetric,
    is_symmetric, compute_violation_cost
)


class TestHelperFunctions(unittest.TestCase):
    def test_make_symmetric(self):
        m = 11
        S = {1, 3}
        sym = make_symmetric(S, m)
        self.assertEqual(sym, {1, 3, 8, 10})  # -1=10, -3=8 mod 11

    def test_is_symmetric(self):
        m = 11
        self.assertTrue(is_symmetric({1, 10, 3, 8}, m))
        self.assertFalse(is_symmetric({1, 3}, m))

    def test_Delta(self):
        A = {0, 1, 2}
        B = {0, 1, 2}
        m = 5
        # Delta(A, B, d=1, m=5): count a in A where (a-1) mod 5 in B
        # a=0: -1 mod 5 = 4, not in B
        # a=1: 0 in B, yes
        # a=2: 1 in B, yes
        self.assertEqual(Delta(A, B, 1, m), 2)

    def test_Sigma(self):
        A = {0, 1, 2}
        B = {0, 1, 2}
        m = 5
        # Sigma(A, B, d=2, m=5): count a in A where (2-a) mod 5 in B
        # a=0: 2 in B, yes
        # a=1: 1 in B, yes
        # a=2: 0 in B, yes
        self.assertEqual(Sigma(A, B, 2, m), 3)


class TestBlockCirculantGraph(unittest.TestCase):
    def test_init_enforces_symmetry(self):
        G = BlockCirculantGraph(n=6, D11={1}, D12={0, 1}, D22={2})
        m = G.m  # 11
        self.assertIn(1, G.D11)
        self.assertIn(m - 1, G.D11)  # -1 mod 11 = 10
        self.assertIn(2, G.D22)
        self.assertIn(m - 2, G.D22)  # -2 mod 11 = 9

    def test_adjacency_V1(self):
        G = BlockCirculantGraph(n=6, D11={1, 2}, D12=set(), D22=set())
        # m = 11, so D11 = {1, 2, 9, 10}
        self.assertTrue(G.adjacent(0, 1))   # diff = 1
        self.assertTrue(G.adjacent(0, 2))   # diff = 2
        self.assertFalse(G.adjacent(0, 3))  # diff = 3
        self.assertTrue(G.adjacent(0, 9))   # diff = 9 (= -2 mod 11)
        self.assertTrue(G.adjacent(0, 10))  # diff = 10 (= -1 mod 11)

    def test_adjacency_cross_block(self):
        G = BlockCirculantGraph(n=6, D11=set(), D12={0, 5}, D22=set())
        # m = 11, V1 = {0..10}, V2 = {11..21}
        # Edge (u, v) with u in V1, v in V2: (v - u) mod 11 in D12
        self.assertTrue(G.adjacent(0, 11))   # diff = 11 mod 11 = 0
        self.assertTrue(G.adjacent(0, 16))   # diff = 16 mod 11 = 5
        self.assertFalse(G.adjacent(0, 12))  # diff = 12 mod 11 = 1

    def test_adjacency_matrix(self):
        G = BlockCirculantGraph(n=4, D11={1}, D12={0}, D22={1})
        adj = G.to_adjacency_matrix()
        # Check symmetry
        for u in range(G.N):
            for v in range(G.N):
                self.assertEqual(adj[u][v], adj[v][u])
        # Verify a few edges
        self.assertEqual(adj[0][1], 1)  # D11 contains 1
        self.assertEqual(adj[0][2], 0)  # D11 doesn't contain 2


class TestCommonNeighbors(unittest.TestCase):
    def setUp(self):
        # Create a test graph with n=6
        random.seed(42)
        self.n = 6
        m = 2 * self.n - 1  # 11

        # Create random symmetric sets
        D11_half = set(random.sample(range(1, (m + 1) // 2 + 1), 3))
        D22_half = set(random.sample(range(1, (m + 1) // 2 + 1), 3))
        D12 = set(random.sample(range(m), 5))

        self.G = BlockCirculantGraph(n=self.n, D11=D11_half, D12=D12, D22=D22_half)

    def test_common_neighbors_V1V1_formula_vs_brute(self):
        """Compare formula-based count with brute force for V1-V1 edges."""
        G = self.G
        m = G.m

        for d in range(1, m):
            formula_count = common_neighbors_V1V1(G, d)
            # Pick representative edge: (0, d) in V1
            brute_count = common_neighbors_brute_force(G, 0, d)
            self.assertEqual(formula_count, brute_count,
                            f"V1V1 mismatch at d={d}: formula={formula_count}, brute={brute_count}")

    def test_common_neighbors_V2V2_formula_vs_brute(self):
        """Compare formula-based count with brute force for V2-V2 edges."""
        G = self.G
        m = G.m

        for d in range(1, m):
            formula_count = common_neighbors_V2V2(G, d)
            # Pick representative edge: (m, m+d) in V2
            u, v = m, m + d
            brute_count = common_neighbors_brute_force(G, u, v)
            self.assertEqual(formula_count, brute_count,
                            f"V2V2 mismatch at d={d}: formula={formula_count}, brute={brute_count}")

    def test_common_neighbors_V1V2_formula_vs_brute(self):
        """Compare formula-based count with brute force for V1-V2 edges."""
        G = self.G
        m = G.m

        for d in range(m):
            formula_count = common_neighbors_V1V2(G, d)
            # Pick representative edge: (0, m+d) where 0 in V1, m+d in V2
            u, v = 0, m + d
            brute_count = common_neighbors_brute_force(G, u, v)
            self.assertEqual(formula_count, brute_count,
                            f"V1V2 mismatch at d={d}: formula={formula_count}, brute={brute_count}")

    def test_common_neighbors_multiple_random_graphs(self):
        """Test with multiple random graphs to ensure robustness."""
        for seed in range(10):
            random.seed(seed)
            n = random.randint(4, 8)
            m = 2 * n - 1

            D11_half = set(random.sample(range(1, (m + 1) // 2 + 1),
                                         random.randint(1, (m + 1) // 4)))
            D22_half = set(random.sample(range(1, (m + 1) // 2 + 1),
                                         random.randint(1, (m + 1) // 4)))
            D12 = set(random.sample(range(m), random.randint(1, m // 2)))

            G = BlockCirculantGraph(n=n, D11=D11_half, D12=D12, D22=D22_half)

            # Test a few random edges
            for _ in range(5):
                d = random.randint(1, m - 1)
                self.assertEqual(common_neighbors_V1V1(G, d),
                               common_neighbors_brute_force(G, 0, d))
                self.assertEqual(common_neighbors_V2V2(G, d),
                               common_neighbors_brute_force(G, m, m + d))

            for _ in range(5):
                d = random.randint(0, m - 1)
                self.assertEqual(common_neighbors_V1V2(G, d),
                               common_neighbors_brute_force(G, 0, m + d))


class TestVerification(unittest.TestCase):
    def test_random_graph_usually_fails(self):
        """Random graphs should almost always violate constraints."""
        random.seed(123)
        num_valid = 0
        num_tests = 20

        for _ in range(num_tests):
            n = 6
            m = 2 * n - 1

            D11_half = set(random.sample(range(1, (m + 1) // 2 + 1),
                                         random.randint(1, (m + 1) // 4)))
            D22_half = set(random.sample(range(1, (m + 1) // 2 + 1),
                                         random.randint(1, (m + 1) // 4)))
            D12 = set(random.sample(range(m), random.randint(1, m // 2)))

            G = BlockCirculantGraph(n=n, D11=D11_half, D12=D12, D22=D22_half)
            result = verify_construction(G)
            if result.valid:
                num_valid += 1

        # Expect very few (if any) random graphs to be valid
        self.assertLess(num_valid, num_tests // 2,
                       f"Too many random graphs are valid: {num_valid}/{num_tests}")

    def test_violation_cost_zero_iff_valid(self):
        """Violation cost should be 0 exactly when construction is valid."""
        random.seed(456)

        for _ in range(10):
            n = 5
            m = 2 * n - 1
            D11 = set(random.sample(range(1, m), random.randint(1, 3)))
            D12 = set(random.sample(range(m), random.randint(1, 3)))
            D22 = set(random.sample(range(1, m), random.randint(1, 3)))

            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
            result = verify_construction(G)
            cost = compute_violation_cost(G)

            if result.valid:
                self.assertEqual(cost, 0)
            else:
                self.assertGreater(cost, 0)


class TestBlueNeighborFormula(unittest.TestCase):
    """Test the blue common neighbor formula (inclusion-exclusion)."""

    def blue_common_brute(self, G, u, v):
        """Count vertices that are NOT adjacent to either u or v (brute force)."""
        count = 0
        for w in range(G.N):
            if w != u and w != v:
                if not G.adjacent(u, w) and not G.adjacent(v, w):
                    count += 1
        return count

    def test_blue_formula_matches_brute_force(self):
        """Verify blue common neighbor formula matches brute force counting."""
        random.seed(789)

        for seed in range(10):
            random.seed(seed)
            n = random.randint(4, 8)
            m = 2 * n - 1
            N = 4 * n - 2

            D11_half = set(random.sample(range(1, (m + 1) // 2 + 1),
                                         random.randint(1, (m + 1) // 4)))
            D22_half = set(random.sample(range(1, (m + 1) // 2 + 1),
                                         random.randint(1, (m + 1) // 4)))
            D12 = set(random.sample(range(m), random.randint(1, m // 2)))

            G = BlockCirculantGraph(n=n, D11=D11_half, D12=D12, D22=D22_half)

            d1 = len(G.D11) + len(G.D12)
            d2 = len(G.D22) + len(G.D12)

            # Test V1-V1 blue edges
            for d in range(1, m):
                if d not in G.D11:
                    red_common = common_neighbors_V1V1(G, d)
                    blue_formula = (N - 2) - d1 - d1 + red_common
                    blue_brute = self.blue_common_brute(G, 0, d)
                    self.assertEqual(blue_formula, blue_brute,
                                    f"V1V1 blue mismatch at d={d}, n={n}: "
                                    f"formula={blue_formula}, brute={blue_brute}")

            # Test V2-V2 blue edges
            for d in range(1, m):
                if d not in G.D22:
                    red_common = common_neighbors_V2V2(G, d)
                    blue_formula = (N - 2) - d2 - d2 + red_common
                    blue_brute = self.blue_common_brute(G, m, m + d)
                    self.assertEqual(blue_formula, blue_brute,
                                    f"V2V2 blue mismatch at d={d}, n={n}: "
                                    f"formula={blue_formula}, brute={blue_brute}")

            # Test V1-V2 blue edges
            for d in range(m):
                if d not in G.D12:
                    red_common = common_neighbors_V1V2(G, d)
                    blue_formula = (N - 2) - d1 - d2 + red_common
                    blue_brute = self.blue_common_brute(G, 0, m + d)
                    self.assertEqual(blue_formula, blue_brute,
                                    f"V1V2 blue mismatch at d={d}, n={n}: "
                                    f"formula={blue_formula}, brute={blue_brute}")

    def test_inclusion_exclusion_identity(self):
        """
        Test the identity: for any edge (u,v),
        red_common + blue_common + deg_u + deg_v - 2*adjacent(u,v) = N - 2

        Rearranged: red_common + blue_common = (N - 2) - deg_u - deg_v + 2*adjacent(u,v)
        For blue edges, adjacent(u,v)=0, so: blue_common = (N - 2) - deg_u - deg_v + red_common
        """
        random.seed(999)

        for _ in range(10):
            n = random.randint(4, 7)
            m = 2 * n - 1
            N = 4 * n - 2

            D11 = set(random.sample(range(1, m), random.randint(2, 4)))
            D12 = set(random.sample(range(m), random.randint(2, 4)))
            D22 = set(random.sample(range(1, m), random.randint(2, 4)))

            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)

            # Pick random edges and verify identity
            for _ in range(20):
                u = random.randint(0, N - 1)
                v = random.randint(0, N - 1)
                if u == v:
                    continue

                red_common = common_neighbors_brute_force(G, u, v)
                blue_common = self.blue_common_brute(G, u, v)

                # Count degree of u and v
                deg_u = sum(1 for w in range(N) if w != u and G.adjacent(u, w))
                deg_v = sum(1 for w in range(N) if w != v and G.adjacent(v, w))

                adj_uv = 1 if G.adjacent(u, v) else 0

                # The identity: every vertex w (w != u, w != v) is in exactly one of:
                # - red common (adjacent to both u and v)
                # - blue common (adjacent to neither u nor v)
                # - adjacent to u only
                # - adjacent to v only
                # Total = N - 2
                expected = N - 2
                actual = red_common + blue_common + (deg_u - adj_uv - red_common) + (deg_v - adj_uv - red_common)

                # Simplify: red_common + blue_common + deg_u - adj_uv - red_common + deg_v - adj_uv - red_common
                # = blue_common + deg_u + deg_v - 2*adj_uv - red_common
                # For this to equal N-2: blue_common = (N-2) - deg_u - deg_v + 2*adj_uv + red_common

                self.assertEqual(actual, expected,
                                f"Identity failed for n={n}, u={u}, v={v}")


class TestKnownConstruction(unittest.TestCase):
    """Test against known constructions from the literature."""

    def test_paley_style_n5(self):
        """
        For n=5, m=9. Try a Paley-like construction.
        Quadratic residues mod 9: {1, 4, 7} (1^2=1, 2^2=4, 4^2=7 mod 9)
        Note: 9 is not prime, so this may not work perfectly.
        """
        n = 5
        m = 9
        # Quadratic residues mod 9: {1, 4, 7}
        QR = {1, 4, 7}
        G = BlockCirculantGraph(n=n, D11=QR, D12=QR, D22=QR)
        result = verify_construction(G)

        # Just verify the test runs - may or may not be valid
        self.assertIsNotNone(result)
        print(f"\nPaley-style n=5: valid={result.valid}, "
              f"max_red={result.max_red_common}/{result.red_threshold}, "
              f"max_blue={result.max_blue_common}/{result.blue_threshold}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
