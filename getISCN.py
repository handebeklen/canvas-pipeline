import sys
import argparse

def read_band_file(band_file):
    bands = []
    with open(band_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            chrom, chrom_start, chrom_end, name, gie_stain = line.strip().split()
            bands.append({
                "chrom": chrom,
                "chromStart": int(chrom_start),
                "chromEnd": int(chrom_end),
                "name": name
            })
    return bands

def find_bands(chrom, start, end, bands):
    start_band = None
    end_band = None
    for band in bands:
        if band['chrom'] == chrom:
            if band['chromStart'] <= start <= band['chromEnd']:
                if not start_band or band['chromStart'] < start_band['chromStart']:
                    start_band = band
            if band['chromStart'] <= end <= band['chromEnd']:
                if not end_band or band['chromEnd'] > end_band['chromEnd']:
                    end_band = band
    return start_band, end_band

def process_cnv_file(cnv_file, bands):
    results = []
    with open(cnv_file, 'r') as file:
        for line in file:
            columns = line.strip().split()
            chrom, pos_range = columns[0].split(':')
            start, end = map(int, pos_range.split('-'))
            cn = columns[3].split('=')[-1]  # Get the CN value

            start_band, end_band = find_bands(chrom, start, end, bands)

            if start_band and end_band:
                if start_band['name'] == end_band['name']:
                    iscn = f"{chrom[3:]}{start_band['name']}({start}_{end})x{cn}"
                else:
                    iscn = f"{chrom[3:]}{start_band['name']}{end_band['name']}({start}_{end})x{cn}"
            else:
                iscn = f"{chrom[3:]}({start}_{end})x{cn}"

            results.append(line.strip() + f" {iscn}")
    
    return results

def process_bed_file(bed_file, bands):
    results = []
    with open(bed_file, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) < 4:  # Skip lines that don't have enough columns
                continue
                
            chrom = f"chr{columns[0]}" if not columns[0].startswith('chr') else columns[0]
            start = int(columns[1])
            end = int(columns[2])
            cn = columns[3]  # Assuming CN is in the 4th column
            
            start_band, end_band = find_bands(chrom, start, end, bands)
            
            if start_band and end_band:
                if start_band['name'] == end_band['name']:
                    iscn = f"{chrom[3:]}{start_band['name']}({start}_{end})x{cn}"
                else:
                    iscn = f"{chrom[3:]}{start_band['name']}{end_band['name']}({start}_{end})x{cn}"
            else:
                iscn = f"{chrom[3:]}({start}_{end})x{cn}"
                
            results.append(f"{chrom}:{start}-{end}\tCN={cn}\t{iscn}")
    
    return results

def write_output(results, output_file):
    with open(output_file, 'w') as file:
        file.write("\n".join(results) + "\n")

def main():
    parser = argparse.ArgumentParser(description='Generate ISCN notation for CNV data')
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--cnv', '-c', 
                            help='Input CNV file')
    input_group.add_argument('--bed', '-d',
                            help='Input BED file')
    parser.add_argument('--band', '-b',
                        required=True,
                        help='Cytogenetic band file')
    parser.add_argument('--output', '-o',
                        required=True,
                        help='Output file path')
    
    args = parser.parse_args()

    bands = read_band_file(args.band)
    
    if args.cnv:
        results = process_cnv_file(args.cnv, bands)
    else:
        results = process_bed_file(args.bed, bands)
        
    write_output(results, args.output)
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()

