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
        pfb_df = pd.read_csv(pfb_path, sep=r"\s+", low_memory=False)

        # Create a copy to avoid SettingWithCopyWarning
        manifest_probes = manifest_df[["Name", "Chromosome", "Position"]].copy()
        manifest_probes.columns = ["Name", "Chr", "Position"]

        # Ensure matching data types
        manifest_probes["Chr"] = manifest_probes["Chr"].astype(str)
        manifest_probes["Position"] = pd.to_numeric(
            manifest_probes["Position"], errors="coerce"
        )

        # Find missing probes
        pfb_probes = set(pfb_df["Name"])
        manifest_probes = manifest_probes.dropna()
        manifest_probe_names = set(manifest_probes["Name"])

        # Find probes missing from PFB
        missing_probes = manifest_probes[
            ~manifest_probes["Name"].isin(pfb_probes)
        ].copy()

        # Find probes in PFB but not in manifest
        extra_probes = pfb_df[~pfb_df["Name"].isin(manifest_probe_names)]
        print("\nNumber of probes in PFB but not in manifest:", len(extra_probes))
        print(extra_probes)

        # Add a default frequency value for the missing probes
        missing_probes["PFB"] = 0.001

        # Combine the original PFB data with the missing probes
        updated_pfb_df = pd.concat([pfb_df, missing_probes], ignore_index=True)

        # Sort the updated PFB data by chromosome and position
        updated_pfb_df["Chr"] = pd.to_numeric(updated_pfb_df["Chr"], errors="coerce")
        updated_pfb_df = updated_pfb_df.sort_values(
            by=["Chr", "Position"], ignore_index=True
        )

        # Save the updated PFB file
        updated_pfb_df.to_csv(output_path, sep="\t", index=False, header=False)
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
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-m", "--manifest", required=True, help="Path to manifest CSV file"
    )

    parser.add_argument("-p", "--pfb", required=True, help="Path to input PFB file")

    parser.add_argument(
        "-o", "--output", required=True, help="Path to save updated PFB file"
    )

    args = parser.parse_args()

    process_pfb(args.manifest, args.pfb, args.output)


if __name__ == "__main__":
    main()
