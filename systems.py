import pandas as pd
import numpy as np
from typing import Literal


def appoint_divisor(
    df: pd.DataFrame,
    district_seats: pd.Series,
    type: Literal["dhondt", "sainte-lague"] = "dhondt",
    party_threshold: float = 0.0
) -> pd.DataFrame:
    
    if type == "dhondt":
        denominator = lambda x: x + 1
    elif type == "sainte-lague":
        denominator = lambda x: 2 * x + 1
    else:
        raise ValueError(type)
    
    # copy of the original array, sorting by votes and adding an extra column
    # for results
    votes = df.copy().sort_values(
        by=["district", "pact", "party", "votes"],
        ascending=[True, True, True, False]
    )
    votes["valid"] = True
    votes["elected"] = False
    
    # calculate national vote percentages for parties
    national_votes_party = votes.groupby(
        ["pact", "party"]
    )["votes"].sum()
    percentage_party = national_votes_party / national_votes_party.sum()
    
    # invalidate parties below the threshold (this dataframe will be used
    # for selecting candidates)
    for pact, party in national_votes_party.index:
        
        if percentage_party[(pact, party)] < party_threshold:
            votes.loc[votes[votes["party"] == party].index, "valid"] = False

    for district in district_seats.index:
        
        # seats to allocate in district
        n_seats = district_seats[district]
        
        # candidates in district
        votes_district = votes[votes["district"]==district]
        
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
                chosen_candidate_idx = votes[
                    (votes["district"] == district)
                    & (votes["party"] == chosen_party)
                    & (votes["elected"] == False)
                    & (votes["valid"] == True)
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
            votes.loc[chosen_candidate_idx, "elected"] = True
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
    
    return votes.drop(columns="valid")