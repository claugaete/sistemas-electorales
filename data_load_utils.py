import pandas as pd

# Useful functions for creating and modifying the results dataframe.

def add_percentage(df: pd.DataFrame, column_name="percentage"):
    """
    Adds a column with the percentage of district-level votes obtained by each
    candidate.
    """
    
    district_votes = (
        df.groupby("district")["votes"]
        .sum()
        .rename("district_votes")
    )
    df = df.merge(
        district_votes,
        how="inner",
        left_on="district",
        right_on=district_votes.index
    )
    df[column_name] = (
        df["votes"] / df["district_votes"]
    )
    return df.drop(columns="district_votes")


def redistrict_candidates(df: pd.DataFrame, subdistricts: pd.Series):
    """
    Splits each district of dataframe `df` into the number of subdistricts
    given by `subdistricts`. Each district's candidates are assigned to one of
    the subdistricts by their index value.
    """
    
    df_mod = df.copy()
    df_mod["district"] = ""
    
    max_n_districts = subdistricts.max()
    letters = [chr(i) for i in range(ord("a"), ord("a")+max_n_districts)]
    
    for idx, row in df.iterrows():
        
        suffix = letters[idx % subdistricts[row["district"]]]
        df_mod.loc[idx, "district"] = str(row["district"]) + suffix
    
    return add_percentage(df_mod)