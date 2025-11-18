import os
import sys
import json

def nc_dataset_map():
    nc_taxid_name_map = {
                            "dengue/all": "3052464",
                            "flu/h1n1pdm": "2955291",
                            "flu/h3n2": "2955291",
                            "flu/vic": "2955465",
                            "sars-cov-2": "3418604",
                            "herpes": "3050292",
                            "hmpv": "3048148",
                            "measles": "3052345",
                            "mpox": "3431483",
                            "mumps": "3052560",
                            "orthoebolavirus": "3052462",
                            "rsv": "3049954",
                            "rubella": "2846071",
                            "yellow-fever": "3046277"
                            }
    return nc_taxid_name_map

def list_leaf_dirs(start):
    leaf_dirs = [p for p, dirs, files in os.walk(start) if not dirs]
    if not leaf_dirs:
        return []
    return leaf_dirs

def construct_nc_index(git_nc_db_dir):
    if not git_nc_db_dir:
        return {}
    nc_taxid_name_map = nc_dataset_map()
    segmented = ["2955291", "2955465"]
    git_data_index = {}
    flu_seg_map = {
        "pb2": "1",
        "pb1": "2",
        "pa": "3",
        "ha": "4",
        "np": "5",
        "na": "6",
        "mp": "7",
        "ns": "8",
    }
    leaf_dirs = list_leaf_dirs(git_nc_db_dir)
    for leaf_dir in leaf_dirs:
        for name, taxid in nc_taxid_name_map.items():
            if f"/{name}" in leaf_dir or f"/{name}/" in leaf_dir:
                if taxid not in segmented:
                    if taxid not in git_data_index:
                        git_data_index[taxid] = {"ALL": [leaf_dir]}
                    else:
                        git_data_index[taxid]["ALL"].append(leaf_dir)
                else:
                    print(f"{leaf_dir = }")
                    for gene, segnum in flu_seg_map.items():
                        print(f"Checking {gene}")
                        if f"/{gene}" in leaf_dir or f"/{gene}/" in leaf_dir:
                            print(f"Found {gene = }")
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
    if not local_nc_db_dir:
        return {}
    refseq_index = {}
    leaf_dirs = list_leaf_dirs(local_nc_db_dir)
    for leaf_dir in leaf_dirs:
        taxid = leaf_dir.split("/")[-2]
        if taxid in refseq_index:
            refseq_index[taxid]["ALL"].append(leaf_dir)
        else:
            refseq_index[taxid] = {"ALL": [leaf_dir]}
    return refseq_index

def write_output(git_nc_db_dir, local_nc_db_dir = None, outfname = "nextclade_index.json", outdir = "."):
    git_data_index = construct_nc_index(git_nc_db_dir) if git_nc_db_dir else {}
    refseq_index = construct_local_index(local_nc_db_dir) if local_nc_db_dir else {}

    combined_index = {}
    combined_index.update(refseq_index)
    combined_index.update(git_data_index)

    abs_outdir = os.path.abspath(outdir)
    outfile = os.path.join(abs_outdir, outfname)

    with open(outfile, "w") as outhandle:
        json.dump(combined_index, outhandle, indent=4, sort_keys = True)

    outhandle.close()

def main():
    git_nc_db_dir = os.path.abspath(sys.argv[1])
    local_nc_db_dir = os.path.abspath(sys.argv[2])
    write_output(git_nc_db_dir, local_nc_db_dir)

if __name__ == "__main__":
    main()
