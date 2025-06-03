#!/usr/bin/env python3
import csv
import sys

def load_ecotypes(ecotypes_file):
    """Load latitude/longitude info from the ecotypes CSV into a dict."""
    latlong = {}
    with open(ecotypes_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Assuming the header has columns: id, latitude, longitude
            latlong[row["id"]] = (row["latitude"], row["longitude"])
    return latlong

def merge_kinship(kinship_file, latlong, output_file):
    """
    Read the kinship file (tab-delimited) and merge lat/long data for both sample pairs.
    Produces an output CSV with additional columns: Lat1, Long1, Lat2, Long2.
    """
    with open(kinship_file, newline='') as kin_f, open(output_file, 'w', newline='') as out_f:
        # Remove any leading '#' from the header line if needed.
        header_line = kin_f.readline().lstrip('#').strip()
        fieldnames = header_line.split()
        # Add new columns for latitude and longitude of both individuals
        fieldnames += ["Lat1", "Long1", "Lat2", "Long2"]
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)
        writer.writeheader()
        
        # Process each line of the kinship file
        for line in kin_f:
            if not line.strip():
                continue
            parts = line.split()
            row = dict(zip(fieldnames[:len(parts)], parts))
            # Lookup lat/long for IID1 and IID2; use "NA" if not found
            iid1 = row.get("IID1")
            iid2 = row.get("IID2")
            lat1, long1 = latlong.get(iid1, ("NA", "NA"))
            lat2, long2 = latlong.get(iid2, ("NA", "NA"))
            row["Lat1"] = lat1
            row["Long1"] = long1
            row["Lat2"] = lat2
            row["Long2"] = long2
            writer.writerow(row)

def main():
    if len(sys.argv) != 4:
        print("Usage: merge_kinship_with_latlong.py <kinship_file> <ecotypes_csv> <output_csv>")
        sys.exit(1)
    kinship_file = sys.argv[1]
    ecotypes_csv = sys.argv[2]
    output_csv = sys.argv[3]
    
    latlong = load_ecotypes(ecotypes_csv)
    merge_kinship(kinship_file, latlong, output_csv)
    print(f"Merged kinship file with lat/long saved to: {output_csv}")

if __name__ == "__main__":
    main()
