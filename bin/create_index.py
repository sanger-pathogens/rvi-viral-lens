import os
import sys
import json

"""
Nextclade datasets Github repo: https://github.com/nextstrain/nextclade_data/tree/master
Nextclade datasets data directory: https://github.com/nextstrain/nextclade_data/tree/master/data

Point this script to the second path in the cloned repo
"""

def nc_dataset_map():
    """
    Map of nextclade dataset key words to corresponding NCBI species taxIDs
    Update this to capture new data
    """
    nc_taxid_name_map = {
                            "dengue/all": "3052464",
                            "flu/h1n1pdm": "2955291",
                            "flu/h3n2": "2955291",
                            "flu/vic": "2955465",
                            "flu/yam": "2955465", ## fluB_yamagata - only HA in nextclade_data
                            "sars-cov-2": "3418604",
                            "herpes": "3050292",
                            "hmpv": "3048148",
                            "measles": "3052345",
                            "mpox": "3431483",
                            "mumps": "3052560",
                            "orthoebolavirus": "3052462",
                            "rsv": "3049954",
                            "rubella": "2846071",
                            "yellow-fever": "3046277",
                            "ev-d68": "3428506"
                            }
    return nc_taxid_name_map

def list_leaf_dirs(start):
    """
    Given a base dir to start at, enumerate paths to deepest (leaf) directories in that start dir
    """
    leaf_dirs = [p for p, dirs, files in os.walk(start) if not dirs]
    if not leaf_dirs:
        return []
    return leaf_dirs

def construct_nc_index(git_nc_db_dir):
    """
    Using nc_dataset_map, enumerate paths to datasets inside the nextclade_data github repo
    """
    nc_taxid_name_map = nc_dataset_map()
    segmented = ["2955291", "2955465"]
    git_data_index = {}
    flu_A_map = {
        "pb2": "1",
        "pb1": "2",
        "pa": "3",
        "ha": "4",
        "np": "5",
        "na": "6",
        "mp": "7",
        "ns": "8",
    }

    flu_B_map = {
        "pb1": "1",
        "pb2": "2",
        "pa": "3",
        "ha": "4",
        "np": "5",
        "na": "6",
        "mp": "7",
        "ns": "8",
    }

    ## collect all leaf directories
    leaf_dirs = list_leaf_dirs(git_nc_db_dir)
    for leaf_dir in leaf_dirs:
        for name, taxid in nc_taxid_name_map.items():
            ## scan for keywords, and if found, add to index dict
            if f"/{name}" in leaf_dir or f"/{name}/" in leaf_dir:
                ## if taxid NOT in list of known segmented viruses, update dict
                if taxid not in segmented:
                    if taxid not in git_data_index:
                        git_data_index[taxid] = {"ALL": [leaf_dir]}
                    else:
                        git_data_index[taxid]["ALL"].append(leaf_dir)
                ## if taxid in list of known segmented viruses
                ## use furter keywords to assigned index entries
                else:
                    flu_seg_map = flu_B_map if "/vic/" in leaf_dir else flu_A_map
                    for gene, segnum in flu_seg_map.items():
                        if f"/{gene}" in leaf_dir or f"/{gene}/" in leaf_dir:
                            if taxid not in git_data_index:
                                git_data_index[taxid] = {segnum: [leaf_dir]}
                            else:
                                if segnum in git_data_index[taxid]:
                                    git_data_index[taxid][segnum].append(leaf_dir)
                                else:
                                    git_data_index[taxid].update({segnum: [leaf_dir]})
                            break
                        continue


    return git_data_index

def construct_local_index(local_nc_db_dir):
    """
    Given the base direectory to the local custom nextclade datasets,
    add them to index.
    NOTE: Does not support segment numbers/subtypes or related directory structures in the current form
    """
    refseq_index = {}
    leaf_dirs = list_leaf_dirs(local_nc_db_dir)
    for leaf_dir in leaf_dirs:
        taxid = leaf_dir.split("/")[-2]
        if taxid in refseq_index:
            refseq_index[taxid]["ALL"].append(leaf_dir)
        else:
            refseq_index[taxid] = {"ALL": [leaf_dir]}
    return refseq_index

def write_output(git_nc_db_dir, local_nc_db_dir = None, nc_include = ["nextstrain"], outfname = "nextclade_index.json", outdir = "."):
    """
    Driver function to create index dictionaries from the local nextclade dir
    and cloned copy of nextclade_data github repo
    """
    combined_index = {}

    ## first create index data from local dir, and add to combined_index
    refseq_index = construct_local_index(local_nc_db_dir)
    combined_index.update(refseq_index)

    ## next, for each subdir in base_path (see below) that the user specified, create an index and update combined_index
    ## base_path refers to: path/to/github/repo/nextclade_data/data
    ## default subdir to include is "nextstrain", ie all data in path/to/github/repo/nextclade_data/data/nextstrain
    for dir_to_include in nc_include:
        git_subdir = os.path.abspath(os.path.join(git_nc_db_dir, dir_to_include))
        if not os.path.exists(git_subdir):
            print(f"Path: {git_subdir} not found. Skipping.")
            continue
        git_data_index = construct_nc_index(git_subdir)
        combined_index.update(git_data_index)

    abs_outdir = os.path.abspath(outdir)
    outfile = os.path.join(abs_outdir, outfname)

    ## dump index to JSON
    with open(outfile, "w") as outhandle:
        json.dump(combined_index, outhandle, indent=4, sort_keys = True)

    outhandle.close()

def main():
    if sys.argv[1] == "-h":
        print("USAGE:")
        print("python create_index.py path/to/nextclade_data/data path/to/local/nextclade/datasets nextstrain,enpen")
        exit(0)

    ## path to cloned version of github repo nextclade_data/data
    ## usually path/to/nextclade_data/data
    git_nc_db_dir = os.path.abspath(sys.argv[1])

    ## path to base dir of local custom nextclade datasets
    local_nc_db_dir = os.path.abspath(sys.argv[2])

    ## subdirs under nextclade_data/data to include
    ## given as a comma-separated list
    nc_include = sys.argv[3].split(",") if sys.argv[3] else ["nextstrain"]
    write_output(git_nc_db_dir, local_nc_db_dir, nc_include)

if __name__ == "__main__":
    main()
