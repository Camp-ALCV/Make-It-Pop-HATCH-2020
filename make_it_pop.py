import pandas as pd
from helper_functions import *
import sys
import os
import argparse
from tqdm.autonotebook import tqdm


def make_bed_df(fullpath):
    print("Reading:", fullpath)
    _f = read_tsv(fullpath)
#     tqdm.pandas(desc="Hash")
#     _f["hash"] = _f.iloc[:,[1,2,4]].progress_apply(lambda x: hash(tuple(x)), axis=1)
    _f["hash"] = _f.iloc[:, 1].astype(str) + _f.iloc[:, 2].astype(str) + _f.iloc[:, 4].astype(str)
#     tqdm.pandas(desc="Chr Hash")
    _f["chr_hash"] = _f.iloc[:, 0].astype(str) + _f.iloc[:, 1].astype(str) + _f.iloc[:, 2].astype(str) + _f.iloc[:, 3].astype(str) + _f.iloc[:, 4].astype(str)
#     _f["chr_hash"] = _f.iloc[:,[0,1,2,3,4]].progress_apply(lambda x: hash(tuple(x)), axis=1)
    return _f

def make_trio_dfs(father: str, mother: str, child: str):
    father = make_bed_df(father)
    mother = make_bed_df(mother)
    child = make_bed_df(child)
    return father, mother, child

def adam_func(father: pd.DataFrame, mother: pd.DataFrame, child: pd.DataFrame):
    exons = pd.read_csv("hg19_exon_locations_with_genes.csv")

    # this might be the more correct way but its slow..
#     set_mother = set(mother["hash"])
#     set_father = set(father["hash"])
#     set_child = set(child["hash"])
#     intercet = set.intersection(set_mother, set_father, set_child)
#     interesting_short = child[child["hash"].isin(intercet)]
    mother_father_hashes = mother[mother["hash"].isin(father["hash"])]["hash"]
    interesting_short = child[~child["hash"].isin(mother_father_hashes)]
#     mother_father_hashes = mother[mother["hash"].isin(father["hash"])]["hash"]
#     interesting_short = son[~son["hash"].isin(mother_father_hashes)]


    interIter = interesting_short.iterrows();
    inter = next(interIter)

    # differences in either parent or just child
    interIndices = []

    for index, row in tqdm(exons.iterrows()):
        while(inter[1].START_POSITION <= row.END_POSITION):
            if (inter[1].START_POSITION >= row.START_POSITION & inter[1].START_POSITION < row.END_POSITION) \
            | \
            (inter[1].START_POSITION >= row.START_POSITION & inter[1].START_POSITION < row.END_POSITION):
                 interIndices.append(inter[0])
            inter = next(interIter)
    return child.iloc[interIndices]

def vicky_func(father: pd.DataFrame, mother: pd.DataFrame, child: pd.DataFrame, lb=0.7, ub=1.04):
    child = child[(child.iloc[:, 5]==0) & (child.iloc[:, 6] > 0)]
    set_child = set(child["chr_hash"])
    set_mother = set(mother["chr_hash"])
    set_father = set(father["chr_hash"])

    family = set.intersection(set_child, set_mother, set_father)
    father = father[father["chr_hash"].isin(family)]
    mother = mother[mother["chr_hash"].isin(family)]
    child = child[child["chr_hash"].isin(family)]

    mother["vratio"] = mother.iloc[:,5] / mother.iloc[:, 6]
    father["vratio"] = father.iloc[:,5] / father.iloc[:, 6]
    return  father[(father["vratio"] > .7) & (father["vratio"] < 1.3)], mother[(mother["vratio"] > lb) & (mother["vratio"] < ub)], father[(father["vratio"] > .7) & (father["vratio"] < 1.3)]

def all_together(father: str, mother: str, child: str, save_path: str):

    father, mother, child = make_trio_dfs(**trio1)

    a_child = adam_func(father, mother, child)

    v_father, v_mother, v_child = vicky_func(father, mother, child)

    vfather = os.path.join(save_path, "v_father.bed")
    vmother = os.path.join(save_path, "v_mother.bed")
    vchild = os.path.join(save_path, "v_child.bed")
    achild = os.path.join(save_path, "a_child.bed")
    v_father.iloc[:, :6].to_csv(vfather, index=False, header=False, sep="\t")
    v_mother.iloc[:, :6].to_csv(vmother, index=False, header=False, sep="\t")
    v_child.iloc[:, :6].to_csv(vchild, index=False, header=False, sep="\t")

    a_child.iloc[:, :6].to_csv(achild, index=False, header=False, sep="\t")
    return vfather, vmother, vchild, achild



if __name__ == "__main__":

    my_parser = argparse.ArgumentParser(description="Convert `*.bed` files to our filtered files")
    my_parser.add_argument("--father", metavar='father',
                           type=str, help="the absolute path to the father's `*.bed` file.")
    my_parser.add_argument("--mother", metavar='mother',
                           type=str, help="the absolute path to the mother's `*.bed` file.")
    my_parser.add_argument("--child", metavar='child',
                           type=str, help="the absolute path to the child's `*.bed` file.")
    args = my_parser.parse_args()
    args = dict(father=args.father, mother=args.mother, child=args.child)
    all_together(save_path=os.getcwd(), **args)

