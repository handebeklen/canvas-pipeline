import argparse

# Define the pseudoautosomal regions (PAR) for the human genome in a dictionary
PAR_REGIONS = {
    "GRCh37": {
        "chrX": [
            (60001, 2699520),
            (154931044, 155260560)
        ],
        "chrY": [
            (10001, 2649520),
            (59034050, 59363566)
        ]
    },
    "GRCh38": {
        "chrX": [
            (10001, 2781479),
            (155701383, 156030895)
        ],
        "chrY": [
            (10001, 2781479),
            (56887903, 57217415)
        ]
    }
}


def is_in_par(chr_num, start, end, genome_version="GRCh37"):
    """Check if a given region falls within the PAR regions based on the genome version."""
    if chr_num in ["chrX", "chrY"]:
        start, end = int(start), int(end)
        par_regions = PAR_REGIONS.get(genome_version, {}).get(chr_num, [])
        for par_start, par_end in par_regions:
            if (par_start <= start <= par_end) or (par_start <= end <= par_end):
                return True
    return False


def transform_data(input_file, output_file, sex_file, genome_version):
    with open(input_file, 'r') as infile, open(sex_file, 'r') as sexfile, open(output_file, 'w') as outfile:
        genders = {line.split()[0]: line.split()[1] for line in sexfile.read().splitlines()}
        num_input = 0
        num_output = 0
        for line in infile:
            num_input += 1
            parts = line.strip().split()
            chr_info, numsnp_info, length_info, state_info, file_info, startsnp_info, endsnp_info, conf = parts
            chr_num, start_end = chr_info.split(':')
            start, end = start_end.split('-')
            cn = int(state_info.split(',')[1].split('=')[1])

            sample_id = file_info.split('.')[0]  # Assuming the sample ID is part of the file name
            sex = genders.get(sample_id, "Unknown")

            if chr_num in ["chrX", "chrY"]:
                if sex == "M":
                    if is_in_par(chr_num, start, end, genome_version):
                        # PAR region behaves like autosomal regions
                        variant_type = "DEL" if cn < 2 else "DUP" if cn > 2 else "UNK"
                    else:
                        # Non-PAR region
                        variant_type = "DEL" if cn < 1 else "DUP" if cn > 1 else "UNK"
                elif sex == "F":
                    variant_type = "DEL" if cn < 2 else "DUP" if cn > 2 else "UNK"
                else:
                    variant_type = "UNK"
            else:
                variant_type = "DEL" if cn < 2 else "DUP" if cn > 2 else "UNK"

            if variant_type != "UNK":
                num_output += 1
                outfile.write(f"{chr_num}\t{start}\t{end}\t{variant_type}\n")
            else:
                print(f"{chr_num}\t{start}\t{end}\t{variant_type}\n")
        print(f"num_input={num_input}, num_output={num_output}")


def main():
    parser = argparse.ArgumentParser(description="Transform genomic data to desired format.")
    parser.add_argument("--input_file", help="The input file containing genomic data.")
    parser.add_argument("--output_file", help="The output file to write transformed data.")
    parser.add_argument("--sex_file", help="Tab-separated genders of the samples")
    parser.add_argument("--genome_version", default="GRCh37", help="Genome version to use for PAR regions (GRCh37 or GRCh38)")

    args = parser.parse_args()

    # Pass the genome_version argument to transform_data
    transform_data(args.input_file, args.output_file, args.sex_file, args.genome_version)
    print(f"Transformed data has been written to {args.output_file}")


if __name__ == "__main__":
    main()
