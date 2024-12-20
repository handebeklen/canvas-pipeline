import argparse
import gzip

def find_closest_probes(chrom, start, end, lrr_data):
    """
    Find the closest probes to CNV boundaries in LRR data
    
    Args:
        chrom (str): Chromosome
        start (int): CNV start position
        end (int): CNV end position
        lrr_data (dict): Dictionary of chromosome positions and their LRR values
    
    Returns:
        tuple: (left_probe, right_probe) positions
    """
    if chrom not in lrr_data:
        return None, None
    
    positions = sorted(lrr_data[chrom])
    
    # Find closest left probe
    left_probe = None
    if positions:
        if start <= positions[0]:
            left_probe = positions[0]
        else:
            for pos in positions:
                if pos <= start:
                    left_probe = pos
                    break
    
    # Find closest right probe
    right_probe = None
    if positions:
        if end >= positions[-1]:
            right_probe = positions[-1]
        else:
            for pos in reversed(positions):
                if pos >= end:
                    right_probe = pos
                    break
    
    return left_probe, right_probe

def read_lrr_file(lrr_file):
    """
    Read gzipped LRR file and organize data by chromosome
    
    Returns:
        dict: Dictionary with chromosomes as keys and sets of positions as values
    """
    lrr_data = {}
    with gzip.open(lrr_file, 'rt') as f:
        for line in f:
            chrom, start, _, _ = line.strip().split()
            # Convert chromosome number to match BED format if needed
            chrom = f"chr{chrom}" if not chrom.startswith('chr') else chrom
            pos = int(start)
            if chrom not in lrr_data:
                lrr_data[chrom] = set()
            lrr_data[chrom].add(pos)
    return lrr_data

def process_probes(bed_file, lrr_file, output_file):
    """
    Process probe data from bed and LRR files
    
    Args:
        bed_file (str): Path to BED format file containing CNV regions
        lrr_file (str): Path to LRR data file
        output_file (str): Path to output file
    """
    # Read LRR data
    lrr_data = read_lrr_file(lrr_file)
    
    # Process BED file and find closest probes
    with open(bed_file, 'r') as bed_f, open(output_file, 'w') as out_f:
        for line in bed_f:
            chrom, start, end, cn = line.strip().split()
            # Add chr prefix if not present
            chrom = f"chr{chrom}" if not chrom.startswith('chr') else chrom
            start = int(start)
            end = int(end)
            
            left_probe, right_probe = find_closest_probes(chrom, start, end, lrr_data)
            
            if left_probe is not None and right_probe is not None:
                # Write output in BED format with original CNV and probe positions
                out_f.write(f"{chrom}\t{left_probe}\t{right_probe}\t{cn}\n")

def main():
    parser = argparse.ArgumentParser(description='Process probe data for CNV regions')
    
    parser.add_argument('--bed', 
                       required=True,
                       help='Input BED file containing CNV regions')
    
    parser.add_argument('--lrr',
                       required=True,
                       help='Input LRR data file')
    
    parser.add_argument('--output',
                       required=True,
                       help='Output file path')
    
    args = parser.parse_args()
    
    process_probes(args.bed, args.lrr, args.output)

if __name__ == "__main__":
    main()
