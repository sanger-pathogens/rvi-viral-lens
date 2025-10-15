#!/usr/bin/env python3

import csv
import json
import argparse

def safeload_json(json_path):
    with open(json_path, 'r') as f:
        try:
            content = json.load(f)
        except Exception as e:
            print(f"Error reading {json_path}:\n{e}")
    return content

def build_summary_csv_mapping():
    summary_map = {
            "sample_id": "Sample_ID",
            "virus": "Virus_Taxon_ID",
            "report_name": "Virus",
            "virus_name": "Species",
            "selected_taxid": "Reference_Taxon_ID",
            "ref_selected": "Selected_Reference",
            "flu_segment": "Flu_Segment",
            "virus_subtype": "Reference_Subtype",
            "sample_subtype": "Sample_Subtype",
            "percent_positions_exceeding_depth_10": "Percentage_of_Genome_Covered",
            "reads_mapped": "Total_Mapped_Reads",
            "bases_mapped": "Total_Mapped_Bases",
            "longest_non_n_subsequence": "Longest_non_N_segment",
            "percent_non_n_bases": "Percentage_non_N_bases",
            "mutations": "total_mutations",
            "insertions": "n_insertions",
            "deletions": "n_deletions",
            "snps": "n_snps",
            "ti_tv_ratio": "ti_tv_ratio",
            "nc.selected_dataset": "nc.selected_dataset",
            "nc.coverage": "nc.coverage",
            "nc.qc.overallScore": "nc.qc.overallScore",
            "nc.qc.overallStatus": "nc.qc.overallStatus",
            "nc.qc.missingData": "nc.qc.missingData",
            "nc.qc.mixedSites": "nc.qc.mixedSites",
            "nc.qc.privateMutations": "nc.qc.privateMutations",
            "nc.qc.snpClusters": "nc.qc.snpClusters",
            "nc.qc.frameShifts": "nc.qc.frameShifts",
            "nc.qc.stopCodons": "nc.qc.stopCodons",
            "id": "file_prefix"
            }

    return summary_map

def build_nc_mappings():
    nc_mappings = {
                    "nc.selected_dataset": lambda d: d["dataset"],
                    "nc.coverage": lambda d: d["results"][0]["coverage"] or "None",
                    "nc.qc.overallScore": lambda d: d["results"][0]["qc"]["overallScore"],
                    "nc.qc.overallStatus": lambda d: d["results"][0]["qc"]["overallStatus"] or "None",
                    "nc.qc.missingData": lambda d: d["results"][0]["qc"]["missingData"] or "None",
                    "nc.qc.mixedSites": lambda d: d["results"][0]["qc"]["mixedSites"] or "None",
                    "nc.qc.privateMutations": lambda d: d["results"][0]["qc"]["privateMutations"] or "None",
                    "nc.qc.snpClusters": lambda d: d["results"][0]["qc"]["snpClusters"] or "None",
                    "nc.qc.frameShifts": lambda d: d["results"][0]["qc"]["frameShifts"] or "None",
                    "nc.qc.stopCodons": lambda d: d["results"][0]["qc"]["stopCodons"] or "None"
                   }
    return nc_mappings

def write_summary_csv(data_list, output_path):
    fieldnames = ['Sample_ID', 'Virus_Taxon_ID', 'Virus', 'Species', 'Reference_Taxon_ID',
                  'Selected_Reference', 'Flu_Segment', 'Reference_Subtype', 'Sample_Subtype',
                  'Percentage_of_Genome_Covered', 'Total_Mapped_Reads', 'Total_Mapped_Bases',
                  'Longest_non_N_segment', 'Percentage_non_N_bases', 'total_mutations', 'n_insertions',
                  'n_deletions', 'n_snps', 'ti_tv_ratio', 'nc.selected_dataset', 'nc.coverage', 'nc.qc.overallScore',
                  'nc.qc.overallStatus', 'nc.qc.missingData', 'nc.qc.mixedSites', 'nc.qc.privateMutations',
                  'nc.qc.snpClusters', 'nc.qc.frameShifts', 'nc.qc.stopCodons', 'file_prefix']

    ## Write to CSV file
    with open(output_path, mode='w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerows(data_list)

    csvfile.close()

def main():
    parser = argparse.ArgumentParser(description="Sort JSON files by primary and fallback key.")
    parser.add_argument("-i", "--sample_id", required=False, default = None, help="Sample ID.")
    parser.add_argument("-m", "--meta_str", required=False, default = None, help="Meta stringified for this script.")
    parser.add_argument("-p", "--qc_json", required=False, default = None, type=str, help="QC properties JSON file path")
    parser.add_argument("-n", "--nc_json", required=False, default = None, type=str, help="Nextclade properties JSON file path")
    parser.add_argument("-rm", "--run_mode", required=False, default="sample", help="Whether the inputs are per-sample or per-run. Allowed values: ['sample', 'run']")
    parser.add_argument("-c", "--concat_files", nargs="+", required=False, default = [], help="List of per-consensus properties files")
    parser.add_argument("-r", "--max_records", required=False, default=1, type=int, help="Number of nextclade datasets to retain in the run-level file.")

    args = parser.parse_args()
    if args.run_mode == "sample":
        if any(i is None for i in [args.sample_id, args.meta_str, args.qc_json]):
            if not args.sample_id:
                print("Missing sample id")
            if not args.meta_str:
                print("Missing meta str")
            if not args.qc_json:
                print("Missing qc json")
            raise ValueError("One or more required argumetns for this mode are missing.")

        meta_dict = safeload_json(args.meta_str)
        qc_map = safeload_json(args.qc_json)
        if args.nc_json:
            nc_list = safeload_json(args.nc_json)
        else:
            nc_list = {}

        sample_merged_json_fname = f"{args.sample_id}.properties.json"

        meta_w_qc = {}
        meta_w_qc.update(meta_dict)
        meta_w_qc.update(qc_map)

        nc_data_mappings = build_nc_mappings()
        if nc_list:
            for k, v_func in nc_data_mappings.items():
                try:
                    val = v_func(nc_list[0])
                    if val == None:
                        val = "None"
                except IndexError as ie:
                    val = "NoResults"
                meta_w_qc.update({k: val})

            meta_w_qc["num_nextclade_datasets"] = len(nc_list)
            meta_w_qc["nextclade_results"] = nc_list

        else:
            for k in nc_data_mappings.keys():
                meta_w_qc.update({k: "NextcladeNotRun"})

            meta_w_qc["num_nextclade_datasets"] = 0
            meta_w_qc["nextclade_results"] = []

        with open(sample_merged_json_fname, "w") as per_con_prop_file_handle:
            json.dump(meta_w_qc, per_con_prop_file_handle, indent=4)

        per_con_prop_file_handle.close()

    elif args.run_mode == "run":
        if not args.concat_files:
            raise ValueError("No files found to concatenate!")

        run_data = []
        csv_data = []
        for json_file in args.concat_files:
            json_data = safeload_json(json_file)
            exisiting_nc_data = json_data["nextclade_results"]
            to_keep = exisiting_nc_data[:args.max_records]
            json_data["nextclade_results"] = to_keep
            csv_data.append(json_data)
            json_data.pop("num_nextclade_datasets", None)
            run_data.append(json_data)

        print(f"{len(run_data)} records in list run_data")
        sorted_run_data_list = sorted(run_data, key= lambda d: (d["sample_id"], d["virus"], d["flu_segment"], d["selected_taxid"], d["nc.qc.overallScore"]))

        run_level_json_fname = "consensus_sequence_properties.json"
        with open(run_level_json_fname, "w") as run_json_handle:
            json.dump(sorted_run_data_list, run_json_handle, indent=4)


        summary_map = build_summary_csv_mapping()
        csv_lines = []
        for per_con_data in csv_data:
            summary_line = {}
            for k, v in per_con_data.items():
                try:
                    summary_line[summary_map[k]] = v
                except KeyError as ke:
                    continue
                    print(f"Key {updated_blob_k} not in summary map, will not be used as a summary report column.")
            csv_lines.append(summary_line)

        ## probably redundant
        sorted_csv_lines = sorted(csv_lines, key=lambda d: (d["Sample_ID"], d["Virus_Taxon_ID"], d["Flu_Segment"], d["Reference_Taxon_ID"], d["nc.qc.overallScore"]))

        summary_csv_fname = "summary_report.csv"

        write_summary_csv(sorted_csv_lines, summary_csv_fname)

if __name__ == "__main__":
    main()
