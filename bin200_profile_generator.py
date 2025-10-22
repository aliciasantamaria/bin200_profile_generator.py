#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import pybedtools
from tqdm import tqdm
from collections import defaultdict

def generate_200bp_binary(peaks_file, input_dir_F, input_dir_M, output_dir_F, output_dir_M):
    # Read the peaks file and validate its format
    peaks = pd.read_csv(
        peaks_file,
        sep='\t',
        skiprows=1,
        header=None,
        names=['chr', 'start', 'end'],
        dtype={'chr': str, 'start': int, 'end': int}
    )
    
    # Basic input checks
    if peaks.empty:
        raise ValueError("Empty peaks file")
    if (peaks['start'] > peaks['end']).any():
        raise ValueError("Invalid peaks with start > end")
    
    chrom = peaks['chr'].unique()[0]
    max_end = peaks['end'].max()

    # Generate 200 bp genomic windows covering the full chromosome
    windows = []
    current_start = 0
    while current_start <= max_end:
        current_end = min(current_start + 199, max_end)
        windows.append((chrom, current_start, current_end))
        current_start += 200

    # Convert windows to DataFrame and index them
    windows_df = pd.DataFrame(windows, columns=['chr', 'start', 'end'])
    windows_dict = {(row['chr'], row['start'], row['end']): idx 
                    for idx, row in windows_df.iterrows()}

    # Convert to BedTool objects for overlap operations
    peaks_bed = pybedtools.BedTool.from_dataframe(peaks)
    windows_bed = pybedtools.BedTool.from_dataframe(windows_df[['chr', 'start', 'end']])

    # Inner function to process each sample directory
    def process_samples(input_dir, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        sample_files = [f for f in os.listdir(input_dir) 
                       if f.endswith(f"_chr{chrom}_binary.txt")]
        
        # Iterate through all binary files for this chromosome
        for file in tqdm(sample_files, desc=f"Processing {os.path.basename(input_dir)}"):
            input_path = os.path.join(input_dir, file)
            output_path = os.path.join(output_dir, file)
            
            # Read header and binary data
            with open(input_path) as f:
                header = [next(f) for _ in range(2)]
                binaries = [int(line.strip()) for line in f]
            
            # Skip if number of binary entries doesn’t match number of peaks
            if len(binaries) != len(peaks):
                print(f"Skipping {file} (binary count mismatch: {len(binaries)} vs {len(peaks)})")
                continue
            
            # Attach binary values to peaks
            peaks_with_binary = peaks.copy()
            peaks_with_binary['binary'] = binaries
            
            # Find overlaps between windows and peaks
            overlaps = windows_bed.intersect(
                pybedtools.BedTool.from_dataframe(peaks_with_binary), 
                wa=True, 
                wb=True, 
                loj=True
            )
            
            # Initialize result list with 0s, one per window
            results = [0] * len(windows_df)
            
            # For each overlap, assign the binary value to the corresponding window
            for interval in overlaps:
                a_chr = interval.fields[0]
                a_start = int(interval.fields[1])
                a_end = int(interval.fields[2])
                window_idx = windows_dict[(a_chr, a_start, a_end)]
                
                # If there is a valid overlap, update binary value
                if len(interval.fields) >= 7 and interval.fields[6] != '.':
                    binary = int(interval.fields[6])
                    results[window_idx] = max(results[window_idx], binary)

            # Write only binary values (plus original header)
            with open(output_path, 'w') as f:
                f.writelines(header)
                for val in results:
                    f.write(f"{val}\n")

    # Debug info: list all input files detected
    sample_files = [f for f in os.listdir(input_dir_F) if f.endswith(f"_chr{chrom}_binary.txt")]
    print(f"Files found in {input_dir_F}: {sample_files}") 

    # Process both female and male datasets
    process_samples(input_dir_F, output_dir_F)
    process_samples(input_dir_M, output_dir_M)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 200 bp binary genomic files')
    parser.add_argument('peaks_file', help='BED file with genomic peaks')
    parser.add_argument('--input_dir_F', required=True, help='Directory with female binary input files')
    parser.add_argument('--input_dir_M', required=True, help='Directory with male binary input files')
    parser.add_argument('--output_dir_F', required=True, help='Output directory for female 200 bp profiles')
    parser.add_argument('--output_dir_M', required=True, help='Output directory for male 200 bp profiles')

    args = parser.parse_args()
    
    # Expand home directories
    args.input_dir_F = os.path.expanduser(args.input_dir_F)
    args.input_dir_M = os.path.expanduser(args.input_dir_M)
    args.output_dir_F = os.path.expanduser(args.output_dir_F)
    args.output_dir_M = os.path.expanduser(args.output_dir_M)
    
    # Run the main function
    generate_200bp_binary(
        args.peaks_file,
        args.input_dir_F,
        args.input_dir_M,
        args.output_dir_F,
        args.output_dir_M
    )
