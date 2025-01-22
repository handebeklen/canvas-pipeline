import pandas as pd
import argparse
import sys

def process_pfb(manifest_path, pfb_path, output_path):
    """
    Process PFB file and update it with manifest information.
    
    Args:
        manifest_path (str): Path to manifest CSV file
        pfb_path (str): Path to input PFB file
        output_path (str): Path to save updated PFB file
    """
    try:
        # Load the files
        manifest_df = pd.read_csv(manifest_path, low_memory=False)
        pfb_df = pd.read_csv(pfb_path, sep=r'\s+', low_memory=False)

        # Print initial info
        print(f"\nTotal probes in manifest: {len(manifest_df)}")
        print(f"Total probes in PFB: {len(pfb_df)}")

        # Create a mapping dictionary from manifest for chromosome and position
        manifest_info = manifest_df.set_index('Name')[['Chromosome', 'Position']].to_dict('index')

        # Update PFB data with manifest information
        for idx, row in pfb_df.iterrows():
            probe_name = row['Name']
            if probe_name in manifest_info:
                pfb_df.at[idx, 'Chr'] = manifest_info[probe_name]['Chromosome']
                pfb_df.at[idx, 'Position'] = manifest_info[probe_name]['Position']

        # Create missing probes DataFrame
        manifest_probes = manifest_df[["Name", "Chromosome", "Position"]].copy()
        manifest_probes.columns = ["Name", "Chr", "Position"]

        # Find missing probes
        pfb_probes = set(pfb_df["Name"])
        manifest_probe_names = set(manifest_probes["Name"])

        # Find probes missing from PFB
        missing_probes = manifest_probes[~manifest_probes["Name"].isin(pfb_probes)].copy()

        # Find probes in PFB but not in manifest
        extra_probes = pfb_df[~pfb_df["Name"].isin(manifest_probe_names)]
        print("\nNumber of probes in PFB but not in manifest:", len(extra_probes))
        print(extra_probes)

        # Add a default frequency value for the missing probes
        missing_probes["PFB"] = 0.001

        # Ensure pfb_df has the correct columns
        pfb_df = pfb_df[["Name", "Chr", "Position", "PFB"]]

        # Combine the original PFB data with the missing probes
        updated_pfb_df = pd.concat([pfb_df, missing_probes], ignore_index=True)

        # Custom sorting function for chromosomes
        def chr_sort(x):
            if pd.isna(x):
                return float('inf')
            try:
                return float(x)
            except ValueError:
                # Handle X, Y, MT etc.
                if x == 'X':
                    return 23
                elif x == 'Y':
                    return 24
                elif x == 'MT':
                    return 25
                else:
                    return float('inf')

        # Sort by chromosome (including X, Y) and position
        updated_pfb_df['Chr_sort'] = updated_pfb_df['Chr'].apply(chr_sort)
        updated_pfb_df = updated_pfb_df.sort_values(by=['Chr_sort', 'Position'], ignore_index=True)
        updated_pfb_df = updated_pfb_df.drop('Chr_sort', axis=1)

        # Ensure final output has exactly these columns in this order
        updated_pfb_df = updated_pfb_df[["Name", "Chr", "Position", "PFB"]]

        # Print final info
        print(f"\nTotal probes in updated PFB: {len(updated_pfb_df)}")
        print(f"Probes with missing chromosome: {updated_pfb_df['Chr'].isna().sum()}")

        # Save the updated PFB file with headers
        updated_pfb_df.to_csv(output_path, sep="\t", index=False)
        print(f"\nUpdated PFB file saved to: {output_path}")
        
    except FileNotFoundError as e:
        print(f"Error: Could not find file - {e}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("Error: One of the input files is empty", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: An unexpected error occurred - {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Update PFB file with manifest information and identify missing/extra probes.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "-m", "--manifest",
        required=True,
        help="Path to manifest CSV file"
    )
    
    parser.add_argument(
        "-p", "--pfb",
        required=True,
        help="Path to input PFB file"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to save updated PFB file"
    )

    args = parser.parse_args()
    
    process_pfb(args.manifest, args.pfb, args.output)

if __name__ == "__main__":
    main()