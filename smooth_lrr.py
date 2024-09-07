#!/usr/bin/env python3

import argparse
import csv
from statistics import median
from itertools import islice


def get_smooth_line(input_file, output_file, smooth_num=10):
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        # Read data into a list
        data = [(row[0], int(row[1]), float(row[3])) for row in reader]

        # Group data into chunks of smooth_num
        for i in range(0, len(data), smooth_num):
            chunk = data[i : i + smooth_num]

            # Calculate the median values
            chr_name = chunk[0][0]
            median_start = int(median([row[1] for row in chunk]))
            median_lrr = median([row[2] for row in chunk])

            # Write the result to the output file
            writer.writerow([chr_name, median_start, median_start, median_lrr])


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Process bedGraph data to compute median values for start and lrr."
    )
    parser.add_argument(
        "--input_file", type=str, required=True, help="Input bedGraph file path"
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output file path for smoothed data",
    )
    parser.add_argument(
        "--smooth_num",
        type=int,
        default=10,
        help="Number of rows to group together for smoothing",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    get_smooth_line(
        input_file=args.input_file,
        output_file=args.output_file,
        smooth_num=args.smooth_num,
    )
