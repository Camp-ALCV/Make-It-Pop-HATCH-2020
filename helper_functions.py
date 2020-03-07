import pandas as pd


def read_tsv(path):

    cols =  ["CHROMOSOME_NAME",
          "START_POSITION",
          "END_POSITION",
          "REFERENCE",
          "ALTERNATE",
          "REFERECE_READS",
          "ALTERNATE_READS"]

    df = pd.read_csv(path, delimiter="\t", names=cols)
    return df
