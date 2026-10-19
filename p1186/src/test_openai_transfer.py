"""Controls for the exact periodic formula and the distinction between transfer interfaces."""

from fractions import Fraction
import unittest

from openai_transfer import check, periodic_count, periodic_probe


class TransferTests(unittest.TestCase):
    def test_periodic_transfer_does_not_require_a_safe_coset(self):
        result = check()
        self.assertTrue(result["passed"])
        self.assertFalse(any(row["safe"] for row in result["proper_cosets"]))
        self.assertEqual(result["direct_periodic_transfer"][0]["limiting_density"], "1/16")
        self.assertEqual(result["direct_periodic_transfer"][1]["limiting_density"], "1/66")

    def test_divisibility_boundaries_and_zero_difference(self):
        self.assertEqual(periodic_count(0, 3, 4), 0)
        self.assertEqual(periodic_count(8, 3, 4), 8)
        self.assertEqual(periodic_count(9, 3, 4), 10)
        self.assertEqual(periodic_count(16, 3, 4), 24)
        self.assertEqual(periodic_count(16, 3, 4, include_degenerate=False), 8)
        for n in range(1, 65):
            remainder = n % 8
            self.assertEqual(Fraction(periodic_count(n, 3, 4), n*n),
                             Fraction(1, 16) + Fraction(1, 2*n) + Fraction(remainder*(8-remainder), 16*n*n))

    def test_requires_the_cyclic_hypothesis(self):
        with self.assertRaisesRegex(ValueError, "cyclic AP-free"):
            periodic_probe([0, 0, 0, 1], 3)

    def test_rejects_extra_colors_even_when_they_are_ap_free(self):
        with self.assertRaisesRegex(ValueError, "binary coloring"):
            periodic_probe(list(range(5)), 3)

    def test_rejects_bad_domains(self):
        for args in ((-1, 3, 4), (5, 1, 4), (5, 3, 0), (2.5, 3, 4)):
            with self.assertRaises(ValueError):
                periodic_count(*args)
        with self.assertRaises(ValueError):
            periodic_count(16, 3, 4, include_degenerate=2)

    def test_planted_transfer_defects(self):
        for plant in ("coloring", "interface", "degenerate", "normalization"):
            with self.subTest(plant=plant):
                self.assertFalse(check(plant)["passed"])


if __name__ == "__main__":
    unittest.main()
