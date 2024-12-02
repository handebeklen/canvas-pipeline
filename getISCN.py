import sys

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

def write_output(results, output_file):
    with open(output_file, 'w') as file:
        file.write("\n".join(results) + "\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <cnv_file> <band_file> <output_file>")
        sys.exit(1)

    cnv_file = sys.argv[1]
    band_file = sys.argv[2]
    output_file = sys.argv[3]

    bands = read_band_file(band_file)
    results = process_cnv_file(cnv_file, bands)
    write_output(results, output_file)
    print(f"Output written to {cnv_file.replace('.txt', '_iscn.txt')}")

if __name__ == "__main__":
    main()

