import pandas as pd
import numpy as np
from typing import Literal


def appoint_divisor(
    df: pd.DataFrame,
    district_seats: pd.Series,
    type: Literal["dhondt", "sainte-lague"] = "dhondt",
    party_threshold: float = 0.0,
    selection_criteria: tuple[str, bool] = ("votes", False),
    reset_votes: bool = True
) -> pd.DataFrame:
    """
    Receives:
    - `df`: a dataframe with candidate information; it must have at least 4
        columns: `district`, `pact`, `party`, and `votes`.
    - `district_seats`: a series containing the number of new seats to be
        assigned to each district.
    - `type`: whether the assignment should be done using D'Hondt/Jefferson or
        Sainte-Laguë/Webster.
    - `party_threshold`: minimum percentage of votes a party must receive to be
        eligible for seat allocation. Must be between 0 and 1.
    - `selection_criteria`: criteria by which the candidates for each party
        will be sorted and selected. First element of the tuple is the column
        name, second element is whether the sort is ascending or not. Defaults
        to sorting by the number of votes (descending).
    - `reset_votes`: whether to create an `elected` column (or modify the
        existing one) so that every candidate starts as not elected. If this
        argument is False, then the `elected` column must already exist in the
        dataframe, and already elected candidates will be ignored.
    
    Returns: a copy of `df`, where the `elected` column says whether each
    candidate was elected to the parliament or not.
    """
    
    if type == "dhondt":
        denominator = lambda x: x + 1
    elif type == "sainte-lague":
        denominator = lambda x: 2 * x + 1
    else:
        raise ValueError(type)
    
    # copy of the original array, sorting by votes and adding an extra column
    # for results
    results = df.copy().sort_values(
        by=["district", "pact", "party", selection_criteria[0]],
        ascending=[True, True, True, selection_criteria[1]]
    )
    results["valid"] = True
    if reset_votes:
        results["elected"] = False
    
    # calculate national vote percentages for parties
    national_votes_party = results.groupby(
        ["pact", "party"]
    )["votes"].sum()
    percentage_party = national_votes_party / national_votes_party.sum()
    
    # invalidate parties below the threshold (this dataframe will be used
    # for selecting candidates)
    for pact, party in national_votes_party.index:
        
        if percentage_party[(pact, party)] < party_threshold:
            results.loc[
                results[results["party"] == party].index, "valid"
            ] = False

    for district in district_seats.index:
        
        # seats to allocate in district
        n_seats = district_seats[district]
        
        # candidates in district
        votes_district = results[results["district"]==district]
        
        # votes for each pact and party
        votes_pact = votes_district.groupby("pact")["votes"].sum()
        votes_party = votes_district.groupby(["pact", "party"])["votes"].sum()
        
        # d'hondt coefficient for each pact and party
        coef_pact = votes_pact.copy().astype(float)
        coef_party = votes_party.copy().astype(float)
        
        # number of seats for each pact and party
        elected_pact = pd.Series([0]*len(votes_pact), index=votes_pact.index)
        elected_party = pd.Series(
            [0]*len(votes_party),
            index=votes_party.index.levels[1]
        )
        
        while n_seats > 0:
            
            # choose the pact with highest coefficient
            chosen_pact_idx = coef_pact.argmax()
            chosen_pact = coef_pact.index[chosen_pact_idx]
            
            # choose the party in the pact with highest coefficient
            chosen_party_idx = coef_party[chosen_pact].argmax()
            chosen_party = coef_party[chosen_pact].index[chosen_party_idx]
            
            try:
                # choose the candidate with the most votes in the party
                # (among the candidates whose party is above the threshold)
                chosen_candidate_idx = results[
                    (results["district"] == district)
                    & (results["party"] == chosen_party)
                    & (results["elected"] == False)
                    & (results["valid"] == True)
                ].iloc[0, :].name
            except:
                # if there are valid candidates left in the party, set its
                # coefficient to 0
                coef_party[(chosen_pact, chosen_party)] = 0
                
                # if all pacts in a party have its coefficient at 0, then set
                # the pact coefficient at 0 as well (they have run out of valid
                # candidates)
                if np.all(coef_party[chosen_pact] == 0):
                    coef_pact[chosen_pact] = 0
            
                # if all pacts have their coefficients at 0, then we have
                # ran out of candidates to pick! simply break
                if np.all(coef_pact == 0):
                    break
                
                continue

            
            # add the seat
            results.loc[chosen_candidate_idx, "elected"] = True
            elected_pact[chosen_pact] += 1
            elected_party[chosen_party] += 1
            
            # update coefficients
            coef_pact[chosen_pact] = (
                votes_pact[chosen_pact]
                / denominator(elected_pact[chosen_pact])
            )
            coef_party[(chosen_pact, chosen_party)] = (
                votes_party[chosen_pact][chosen_party]
                / denominator(elected_party[chosen_party])
            )

            n_seats -= 1
    
    return results.drop(columns="valid")


def appoint_divisor_national(
    df: pd.DataFrame,
    total_seats: int,
    type: Literal["dhondt", "sainte-lague"] = "dhondt",
    party_threshold: float = 0.0,
    selection_criteria: tuple[str, bool] = ("percentage", False),
    reset_votes: bool = True
):
    """
    Appoints seats to pacts and parties on a national level using some divisor
    method (D'Hondt/Jefferson or Sainte-Laguë/Webster). Seats are assigned
    within each party based on the order given by `selection_criteria`
    (defaults to percentage obtained by each candidate in their district,
    in order to not benefit larger districts).
    
    Receives same parameters as the district-level divisor method, only with
    the total number of seats to allocate instead of the district-level
    distribution.
    """
    
    # save original districts and assign everyone to one "big" district
    candidate_districts = df["district"]
    df_unified = df.copy()
    df_unified["district"] = 1
    
    results = appoint_divisor(
        df_unified,
        pd.Series({1: total_seats}),
        type,
        party_threshold,
        selection_criteria,
        reset_votes
    )
    
    results["district"] = candidate_districts
    
    return results