#!/usr/bin/env python3

import unittest

from mhass_cigar_utils import parse_cigar, walk_cigar


class TestCigarWalker(unittest.TestCase):

    def test_parse_cigar(self):
        self.assertEqual(
            parse_cigar("4=1I3X2D"),
            [
                (4, "="),
                (1, "I"),
                (3, "X"),
                (2, "D"),
            ],
        )

    def test_perfect_match(self):
        result = walk_cigar(
            "10=",
            reference_length=10,
            query_length=10,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["alignment_target_start_min"],
            3,
        )
        self.assertEqual(
            result["alignment_target_start_max"],
            3,
        )

        self.assertEqual(
            result["alignment_target_end_min"],
            7,
        )
        self.assertEqual(
            result["alignment_target_end_max"],
            7,
        )

        self.assertEqual(result["start_shift_min"], 0)
        self.assertEqual(result["end_shift_min"], 0)

        self.assertFalse(
            result["left_boundary_ambiguous"]
        )
        self.assertFalse(
            result["right_boundary_ambiguous"]
        )

    def test_insertion_before_target(self):
        # Reference:
        # 0123456789
        #
        # Target:
        #    [3,7)
        #
        # One query insertion occurs after reference position 1,
        # within the left wrapper.
        result = walk_cigar(
            "1=1I9=",
            reference_length=10,
            query_length=11,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["left_wrapper_insertions"],
            1,
        )

        self.assertEqual(
            result["alignment_target_start_min"],
            4,
        )
        self.assertEqual(
            result["alignment_target_start_max"],
            4,
        )

        self.assertEqual(result["start_shift_min"], 1)
        self.assertEqual(result["end_shift_min"], 0)

    def test_deletion_before_target(self):
        result = walk_cigar(
            "1=1D8=",
            reference_length=10,
            query_length=9,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["left_wrapper_deletions"],
            1,
        )

        self.assertEqual(
            result["alignment_target_start_min"],
            2,
        )
        self.assertEqual(
            result["alignment_target_start_max"],
            2,
        )

        self.assertEqual(result["start_shift_min"], -1)
        self.assertEqual(result["end_shift_min"], 0)

    def test_insertion_exactly_at_left_boundary(self):
        result = walk_cigar(
            "3=1I7=",
            reference_length=10,
            query_length=11,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["left_boundary_insertions"],
            1,
        )

        self.assertEqual(
            result["alignment_target_start_min"],
            3,
        )
        self.assertEqual(
            result["alignment_target_start_max"],
            4,
        )

        self.assertTrue(
            result["left_boundary_ambiguous"]
        )
    
    # test 6
    def test_insertion_exactly_at_right_boundary(self):
        result = walk_cigar(
            "7=1I3=",
            reference_length=10,
            query_length=11,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["right_boundary_insertions"],
            1,
        )

        self.assertEqual(
            result["alignment_target_end_min"],
            7,
        )
        self.assertEqual(
            result["alignment_target_end_max"],
            8,
        )

        self.assertTrue(
            result["right_boundary_ambiguous"]
        )
    
    # test 7
    def test_insertion_inside_right_wrapper(self):
        result = walk_cigar(
            "8=1I2=",
            reference_length=10,
            query_length=11,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["right_wrapper_insertions"],
            1,
        )

        self.assertEqual(
            result["alignment_target_end_min"],
            7,
        )
        self.assertEqual(
            result["alignment_target_end_max"],
            7,
        )

        self.assertEqual(result["end_shift_min"], -1)
        self.assertEqual(result["end_shift_max"], -1)

        self.assertFalse(
            result["right_boundary_ambiguous"]
        )
    
    # test 8
    def test_deletion_inside_right_wrapper(self):
        result = walk_cigar(
            "8=1D1=",
            reference_length=10,
            query_length=9,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["right_wrapper_deletions"],
            1,
        )

        self.assertEqual(
            result["alignment_target_end_min"],
            7,
        )
        self.assertEqual(
            result["alignment_target_end_max"],
            7,
        )

        self.assertEqual(result["end_shift_min"], 1)
        self.assertEqual(result["end_shift_max"], 1)

        self.assertFalse(
            result["right_boundary_ambiguous"]
        )
    
    # test 9
    def test_substitutions_by_region(self):
        result = walk_cigar(
            "1=1X2=1X3=1X1=",
            reference_length=10,
            query_length=10,
            target_start=3,
            target_end=7,
        )

        self.assertEqual(
            result["left_wrapper_substitutions"],
            1,
        )

        self.assertEqual(
            result["target_substitutions"],
            1,
        )

        self.assertEqual(
            result["right_wrapper_substitutions"],
            1,
        )

        self.assertEqual(result["start_shift_min"], 0)
        self.assertEqual(result["end_shift_min"], 0)
    
if __name__ == "__main__":
    unittest.main(verbosity=2)
