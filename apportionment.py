import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from typing import Literal
from parliamentarch_geometry import get_seats_centers


def groupby_sum_and_pivot(df, index, columns, values, column_order = None):
    """
    Group values of a dataframe on two columns (`index` and `columns`), sums
    the values of `values` for each group and pivots it into a dataframe where
    the rows are represented by the values of `index`, and the columns by the
    values of `columns`. New column order can also be changed (or columns can
    be removed) by passing the extra parameter `column_order` with a list of
    the new columns to keep.
    """
    
    pivoted_df = (
        df.groupby([index, columns])[values]
        .sum()
        .to_frame().reset_index()
        .pivot(index=index, columns=columns, values=values)
        .replace(to_replace={np.nan: 0})
    )
    
    if column_order is not None:
        return pivoted_df.reindex(columns=column_order)
    else:
        return pivoted_df


class Apportionment:

    def __init__(
        self,
        results: pd.DataFrame,
        district_fixed_seats: pd.Series,
        colors: pd.Series | None = None,
        party_order: list[str] | None = None
    ):
        """
        Receives:
        - `results`: a dataframe with five columns: `candidate`, `pact`,
            `party`, `district`, and `votes`.
        - `district_fixed_seats`: a series with districts on its indices, 
            and the number of guaranteed seats per district as values.
        - `colors`: a series with parties on its indices, and a
            `matplotlib`-accepted color for each party as values.
        - `party_order`: a list with all the parties in `results`, in any
            desired order. It is used for ordering all frames and graphs with
            parties involved.
        """

        self.results = results.copy()

        # if no party order is given, we default to dataframe order
        real_parties = self.results["party"].drop_duplicates().to_list()
        if party_order is None:
            self.parties = real_parties

        else:
            party_order_set = set(party_order)
            missing_parties = set(real_parties).difference(party_order_set)

            # if there are missing parties in the given order, error
            if len(missing_parties) > 0:
                raise ValueError(
                    "party_order is missing parties present in dataframe:" +
                    str(missing_parties)
                )

            # if there are duplicates, error
            if len(party_order) != len(party_order_set):
                raise ValueError(
                    "party_order contains duplicates"
                )

            # use the given order, removing entries that aren't in the
            # actual dataframe
            self.parties = [
                p for p in party_order if p in real_parties
            ]

        self.districts = (
            self.results["district"]
            .drop_duplicates()
            .sort_values()
            .to_list()
        )

        self.colors = colors

        self.pacts = self.results["pact"].drop_duplicates().to_list()

        # add percentage of votes (of their district) to each candidate
        district_votes = (
            self.results.groupby("district")["votes"]
            .sum()
            .rename("district_votes")
        )
        self.results = self.results.merge(
            district_votes,
            how="inner",
            left_on="district",
            right_on=district_votes.index
        )
        self.results["percentage"] = (
            self.results["votes"] / self.results["district_votes"]
        )
        self.results = self.results.drop(columns="district_votes")

        # district-party data
        self.district_party_votes = groupby_sum_and_pivot(
            self.results, "district", "party", "votes", self.parties
        ).astype(int)
        self.district_party_vote_share = groupby_sum_and_pivot(
            self.results, "district", "party", "percentage", self.parties
        )
        self.district_party_results = groupby_sum_and_pivot(
            self.results, "district", "party", "elected", self.parties
        ).astype(int)

        # district seat data
        self.district_fixed_seats = district_fixed_seats
        self.district_real_seats = (
            self.results.groupby("district")["elected"].sum()
        )

        # party data
        self.party_vote_share = (
            self.district_party_votes.sum(axis=0)
            / self.district_party_votes.sum().sum()
        )
        self.party_results = (
            self.district_party_results
            .sum(axis=0)
        )

        # district-pact data
        self.district_pact_votes = groupby_sum_and_pivot(
            self.results, "district", "pact", "votes"
        ).astype(int)
        self.district_pact_vote_share = groupby_sum_and_pivot(
            self.results, "district", "pact", "percentage"
        )
        self.district_pact_results = groupby_sum_and_pivot(
            self.results, "district", "pact", "elected"
        ).astype(int)
        # pact data
        self.pact_vote_share = (
            self.district_pact_votes.sum(axis=0)
            / self.district_pact_votes.sum().sum()            
        )
        self.pact_results = self.district_pact_results.sum(axis=0)

        self.n_seats = self.party_results.sum()

    # measure disproportionality
    def disproportionality_national(
        self,
        type: Literal["gallagher", "loosemore-hanby", "wasted-pct"],
        use_pacts: bool = False
    ) -> float:
        """
        Returns the desired disproportionality index, calculated at the
        national level.
        """

        if type in ("gallagher", "loosemore-hanby"):
            if use_pacts:
                expected = self.pact_vote_share
                results = self.pact_results.copy()
            else:
                expected = self.party_vote_share
                results = self.party_results.copy()

            # normalizamos los resultados
            results /= results.sum()

            if type == "gallagher":
                return np.sqrt(np.sum(np.square(expected - results)) / 2)

            elif type == "loosemore-hanby":
                return np.sum(np.abs(expected - results)) / 2

        elif type == "wasted-pct":
            # special case, wasted percentage uses district-level votes even
            # when counting at a national scale (it counts the votes which went
            # to parties that did not get any seats in each district)
            if use_pacts:
                votes = self.district_pact_votes
                results = self.district_pact_results
            else:
                votes = self.district_party_votes
                results = self.district_party_results

            wasted_votes = (votes * (results == 0)).sum().sum()
            total_votes = votes.sum().sum()
            return wasted_votes / total_votes

        else:
            raise ValueError(type)

    def disproportionality_districts(
        self,
        type: Literal["gallagher", "loosemore-hanby", "wasted-pct"],
        use_pacts: bool = False
    ) -> pd.Series:
        """
        Returns a series with the desired disproportionality index calculated
        at district level for all districts.
        """

        if use_pacts:
            expected = self.district_pact_vote_share
            results = self.district_pact_results.copy()
        else:
            expected = self.district_party_vote_share
            results = self.district_party_results.copy()

        # normalize results for each district
        results = results.div(results.sum(axis=1), axis=0)

        if type == "gallagher":
            return np.sqrt(np.sum(np.square(expected - results), axis=1) / 2)

        elif type == "loosemore-hanby":
            return np.sum(np.abs(expected - results), axis=1) / 2

        elif type == "wasted-pct":
            return np.sum(expected * (results == 0), axis=1)

        else:
            raise ValueError(type)

    # fragmentation
    def effective_num_parties(self, use_pacts: bool = False) -> float:
        """
        Returns the effective number of parties in parliament, using the
        formula by Laakso and Taagepera.
        """

        if use_pacts:
            return 1 / np.sum(np.square(self.pact_results / self.n_seats))

        return 1 / np.sum(np.square(self.party_results / self.n_seats))

    # quota conditions
    def quota_condition(self, district = None) -> pd.DataFrame:
        """
        Returns a dataframe in which the obtained number of seats for each
        party is compared with the expected number of seats (vote share
        multiplied by total seats), in order to calculate whether each party
        fulfills the quota rule.
        """

        if district is None:
            expected_seats = self.party_vote_share * self.party_results.sum()
            results = self.party_results
        else:
            expected_seats = (
                self.district_party_vote_share.loc[district]
                * self.district_real_seats[district]
            )
            results = self.district_party_results.loc[district]

        return pd.DataFrame({
            "expected": expected_seats,
            "seats": results,
            "meets_quota": (
                (results == np.ceil(expected_seats))
                | (results == np.floor(expected_seats))
            ),
        }, index=self.parties)

    def quota_condition_districts(self) -> pd.DataFrame:
        """
        Checks the quota rule for each party at district level, and counts
        the number of parties in each district that meet (and don't meet) the
        rule.
        """

        n_parties = len(self.parties)

        n_parties_with_quota = pd.Series(
            [
                self.quota_condition(d)["meets_quota"].sum()
                for d in self.districts
            ],
            index=self.districts
        )

        return pd.DataFrame({
            "parties_with_quota": n_parties_with_quota,
            "parties_without_quota": n_parties - n_parties_with_quota
        }, index=self.districts)

    # plots
    def plot_quota_diff(self, sort_vote_share=False):
        """
        Plots the difference between expected seats and actual seats for each
        party (as calculated by checking the quota rule). Use `sort` to plot
        parties in order of their vote share, overriding the default party
        order.
        """

        quota_df = self.quota_condition()
        
        if sort_vote_share:
            parties = self.party_vote_share.sort_values(
                ascending=False
            ).index.to_list()
            colors = self.colors.reindex(parties)
            quota_df = quota_df.reindex(index=parties)
        else:
            parties = self.parties
            colors = self.colors

        fig, ax = plt.subplots()

        ax.bar(
            parties,
            quota_df["seats"] - quota_df["expected"],
            color=colors,
        )
        ax.set_xticks(parties)
        ax.set_xticklabels(parties, rotation=90)

        ax.set_xlabel("Partido")
        ax.set_ylabel("Dif. entre escaños reales y esperados")
        fig.suptitle("Disproporcionalidad entre votación y escaños")

        return fig

    def plot_parliament(self):
        """
        Plots the parliament results as a hemicycle (à la Wikipedia).
        """

        sorted_seat_centers = sorted(
            get_seats_centers(self.n_seats).items(),
            key=lambda seat: seat[1],
            reverse=True
        )
        seat_centers_array = np.array([
            [coords[0], coords[1]] for coords, _ in sorted_seat_centers
        ])

        color_list = []
        legend_elems = []
        for party, party_seats in self.party_results.items():

            if party_seats == 0:
                continue

            color_list += [self.colors[party]]*party_seats
            legend_elems.append(Line2D(
                [0],
                [0],
                marker='o',
                color="w",
                markerfacecolor=self.colors[party],
                markersize=15,
                label=f"{party} ({party_seats})"
            ))

        fig, ax = plt.subplots()
        fig.set_size_inches(12, 4.8)

        ax.scatter(
            seat_centers_array[:, 0],
            seat_centers_array[:, 1],
            c=color_list,
            # more seats => dots are smaller (up to a point)
            s=np.clip(30_000/(self.n_seats+1), 0, 500),
        )
        ax.text(1, 0.1, self.n_seats, fontsize=50, ha="center")
        ax.set_axis_off()
        ax.set_xlim(0, 2.5)
        ax.legend(handles=legend_elems, title="Partidos")

        return fig

    # summary of statistics
    def summary(self):
        """
        Prints a summary of the most important statistics, and shows useful
        graphs.
        """

        extra_seats = (self.district_real_seats - self.district_fixed_seats)

        extra_seat_distribution = (
            extra_seats
            .groupby(extra_seats)
            .count()
            .to_string()
        )

        gallagher_parties = self.disproportionality_national(
            type="gallagher"
        )*100
        gallagher_pacts = self.disproportionality_national(
            type="gallagher", use_pacts=True
        )*100
        lh_parties = self.disproportionality_national(
            type="loosemore-hanby"
        )*100
        lh_pacts = self.disproportionality_national(
            type="loosemore-hanby", use_pacts=True
        )*100
        wasted_parties = self.disproportionality_national(
            type="wasted-pct"
        )*100
        wasted_pacts = self.disproportionality_national(
            type="wasted-pct", use_pacts=True
        )*100

        district_disp_df = pd.DataFrame({
            'gallagher-partidos': self.disproportionality_districts(
                type="gallagher"
            ),
            'gallagher-pactos': self.disproportionality_districts(
                type="gallagher", use_pacts=True
            ),
            'lh-partidos': self.disproportionality_districts(
                type="loosemore-hanby"
            ),
            'lh-pactos': self.disproportionality_districts(
                type="loosemore-hanby", use_pacts=True
            ),
            'wasted-partidos': self.disproportionality_districts(
                type="wasted-pct"
            ),
            'wasted-pactos': self.disproportionality_districts(
                type="wasted-pct", use_pacts=True
            )
        }).describe().loc[["mean", "std", "min", "50%", "max"]]

        parties_meet_quota = self.quota_condition()["meets_quota"].sum()
        district_parties_with_quota = self.quota_condition_districts()[
            "parties_with_quota"
        ]
        districts_meet_quota = (
            district_parties_with_quota == len(self.parties)
        ).sum()

        print(
            f"""RESUMEN DE ESTADÍSTICAS:

Número de distritos por cantidad de escaños extra:
N° escaños extra / N° distritos
{extra_seat_distribution}

Disproporcionalidad nacional:
- Gallagher (partidos): {gallagher_parties:.2f}%
- Gallagher (pactos): {gallagher_pacts:.2f}%
- Loosemore-Hanby (partidos): {lh_parties:.2f}%
- Loosemore-Hanby (pactos): {lh_pacts:.2f}%
- Wasted vote percentage (partidos): {wasted_parties:.2f}%
- Wasted vote percentage (pactos): {wasted_pacts:.2f}%

Disproporcionalidad por distrito:
{district_disp_df}

Fragmentación:
- Número de partidos: {(self.party_results > 0).sum()}
- Número efectivo de partidos: {self.effective_num_parties():.2f}
- Número de pactos: {(self.pact_results > 0).sum()}
- Número efectivo de pactos: {self.effective_num_parties(use_pacts=True):.2f}

Quotas:
Partidos que cumplen quota: {parties_meet_quota}/{len(self.parties)}
Distritos donde todos los partidos cumplen quota: {
    districts_meet_quota
}/{len(self.districts)}
"""
        )

        fig_quota = self.plot_quota_diff(sort_vote_share=True)
        fig_parliament = self.plot_parliament()
        plt.show()
