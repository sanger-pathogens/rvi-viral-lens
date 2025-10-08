#!/usr/bin/env python3

import os
import glob
import json
import argparse

def load_json_files(file_paths_list):
    """Given a list of JSON filepaths, read them into a list
    """
    data = []
    for path in file_paths_list:
        assembly_name = path.split("/")[-1].split(".")[1]
        with open(path, 'r') as f:
            try:
                content = json.load(f)
                content["dataset"] = assembly_name
                data.append(content)
            except Exception as e:
                print(f"Error reading {path}:\n{e}")
    return data

def get_primary_key(obj):
    """Get the primary key from the nextclade JSON
    Here: results[0]["qc"]["overallScore"]
    Note the results[0] works because each JSON is expected to have
    exactly one set of results - corresponding to the consensus fasta
    """
    try:
        return obj["results"][0]["qc"]["overallScore"]
    except IndexError as ie:
        if len(obj["results"]) == 0:
            seq = obj["errors"][0]["seqName"]
            print(f"For sequence {seq}, found NO RESULTS.")
    return 1.0

def get_secondary_key(obj):
    """Get the secondary key from the nextclade JSON
    Here: results[0]["coverage"]
    Note the results[0] works because each JSON is expected to have
    exactly one set of results - corresponding to the consensus fasta
    """
    try:
        return obj["results"][0]["coverage"]
    except IndexError as ie:
        if len(obj["results"]) == 0:
            seq = obj["errors"][0]["seqName"]
            print(f"For sequence {seq}, found NO RESULTS.")
    return 0.0

def sort_json_list(json_list):
    """Given a list of nexclade outputs (decoded JSON)
    try to sort by overallScore - min to max
    If all overallScore values are identical, use coverage (max to min)
    """
    # Extract primary keys
    primary_keys = [get_primary_key(obj) for obj in json_list]

    # Check if all primary keys are the same
    all_pk_same = all(pk == primary_keys[0] for pk in primary_keys)

    if all_pk_same:
        print("All primary keys are identical — using secondary key: coverage")
        return sorted(json_list, key=get_secondary_key, reverse = True)
    else:
        return sorted(json_list, key=get_primary_key)

def write_output_json(sorted_list, output_path):
    """Dump sorted list of decoded nextclade JSONs to file
    """
    with open(output_path, 'w') as out_handle:
        json.dump(sorted_list, out_handle, indent=4)


def main():
    parser = argparse.ArgumentParser(description="Sort JSON files by primary and fallback key.")
    parser.add_argument("input", nargs='+', help="Input JSON files (can use wildcards).")
    parser.add_argument("-o", "--output", default="sorted_output.json", help="Output JSON file name.")

    args = parser.parse_args()

    # Expand wildcards (e.g., *.json)
    file_paths = []
    for pattern in args.input:
        file_paths.extend(glob.glob(pattern))

    if not file_paths:
        ("No files found!!!.")
        return

    json_data = load_json_files(file_paths)
    sorted_data = sort_json_list(json_data)
    write_output_json(sorted_data, args.output)
    print(f"Sorted JSON written to {args.output}")

if __name__ == "__main__":
    main()
