#!/usr/bin/env python3

import unittest

from mhass_fastq_utils import (
    orient_sequence_and_qualities,
    slice_sequence_and_qualities,
)


class TestMhassFastqUtils(unittest.TestCase):

    def test_as_read_orientation(self):
        sequence = "ACGTTA"
        qualities = [5, 10, 15, 20, 25, 30]

        oriented_sequence, oriented_qualities = (
            orient_sequence_and_qualities(
                sequence,
                qualities,
                "as_read",
            )
        )

        self.assertEqual(
            oriented_sequence,
            "ACGTTA",
        )

        self.assertEqual(
            oriented_qualities,
            [5, 10, 15, 20, 25, 30],
        )

    def test_reverse_complement_orientation(self):
        sequence = "ACGTTA"
        qualities = [5, 10, 15, 20, 25, 30]

        oriented_sequence, oriented_qualities = (
            orient_sequence_and_qualities(
                sequence,
                qualities,
                "revcomp",
            )
        )

        self.assertEqual(
            oriented_sequence,
            "TAACGT",
        )

        self.assertEqual(
            oriented_qualities,
            [30, 25, 20, 15, 10, 5],
        )

    def test_sequence_quality_mismatch_before_orientation(self):
        with self.assertRaises(ValueError):
            orient_sequence_and_qualities(
                "ACGT",
                [10, 20, 30],
                "as_read",
            )

    def test_invalid_orientation(self):
        with self.assertRaises(ValueError):
            orient_sequence_and_qualities(
                "ACGT",
                [10, 20, 30, 40],
                "sideways",
            )

    def test_slice_sequence_and_qualities(self):
        sequence = "AACCGGTT"
        qualities = [
            10, 11, 12, 13,
            14, 15, 16, 17,
        ]

        sliced_sequence, sliced_qualities = (
            slice_sequence_and_qualities(
                sequence,
                qualities,
                2,
                6,
            )
        )

        self.assertEqual(
            sliced_sequence,
            "CCGG",
        )

        self.assertEqual(
            sliced_qualities,
            [12, 13, 14, 15],
        )

    def test_zero_length_slice(self):
        sequence = "ACGT"
        qualities = [10, 20, 30, 40]

        sliced_sequence, sliced_qualities = (
            slice_sequence_and_qualities(
                sequence,
                qualities,
                2,
                2,
            )
        )

        self.assertEqual(sliced_sequence, "")
        self.assertEqual(sliced_qualities, [])

    def test_negative_start_fails(self):
        with self.assertRaises(ValueError):
            slice_sequence_and_qualities(
                "ACGT",
                [10, 20, 30, 40],
                -1,
                3,
            )

    def test_end_beyond_sequence_fails(self):
        with self.assertRaises(ValueError):
            slice_sequence_and_qualities(
                "ACGT",
                [10, 20, 30, 40],
                1,
                5,
            )

    def test_slice_after_reverse_complement(self):
        sequence = "AACCGGTT"
        qualities = [
            10, 11, 12, 13,
            14, 15, 16, 17,
        ]

        oriented_sequence, oriented_qualities = (
            orient_sequence_and_qualities(
                sequence,
                qualities,
                "revcomp",
            )
        )

        sliced_sequence, sliced_qualities = (
            slice_sequence_and_qualities(
                oriented_sequence,
                oriented_qualities,
                2,
                6,
            )
        )

        self.assertEqual(
            oriented_sequence,
            "AACCGGTT",
        )

        self.assertEqual(
            oriented_qualities,
            [17, 16, 15, 14, 13, 12, 11, 10],
        )

        self.assertEqual(
            sliced_sequence,
            "CCGG",
        )

        self.assertEqual(
            sliced_qualities,
            [15, 14, 13, 12],
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
