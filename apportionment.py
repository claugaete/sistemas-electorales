import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from typing import Literal
from parliamentarch_geometry import get_seats_centers


def merge_columns(df, new_names):
    
    return df.rename(columns=new_names).T.groupby(level=0).sum().T

class Apportionment:

    def __init__(
        self,
        district_party_votes: pd.DataFrame,
        district_fixed_seats: pd.Series,
        district_party_results: pd.DataFrame,
        pacts: dict[str, list[str]] | None = None,
        colors: pd.Series | None = None
    ):
        
        self.parties = district_party_results.columns.to_list()
        self.districts = district_party_results.index.to_list()
        self.colors = colors        
        if pacts is None:
            self.pacts = dict(zip(self.parties, self.parties))
        else:
            self.pacts = pacts

        # district-party data
        self.district_party_votes = district_party_votes
        self.district_party_expected = district_party_votes.div(
            district_party_votes.sum(axis=1), axis=0
        )
        self.district_party_results = district_party_results

        # district seat data
        self.district_fixed_seats = district_fixed_seats
        self.district_real_seats = district_party_results.sum(axis=1)

        # party data
        self.party_expected = (
            district_party_votes.sum(axis=0)
            / district_party_votes.sum().sum()
        )
        self.party_results = district_party_results.sum(axis=0)

        party_to_pact = []
        for pact, parties in self.pacts.items():
            for party in parties:
                party_to_pact.append((party, pact))
        # district-pact data
        self.district_pact_votes = merge_columns(
            district_party_votes, dict(party_to_pact)
        )
        self.district_pact_expected = self.district_pact_votes.div(
            self.district_pact_votes.sum(axis=1), axis=0
        )
        self.district_pact_results = merge_columns(
            district_party_results, dict(party_to_pact)
        )
        # pact data
        self.pact_expected = (
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

        if type in ("gallagher", "loosemore-hanby"):
            if use_pacts:
                expected = self.pact_expected
                results = self.pact_results.copy()
            else:
                expected = self.party_expected
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

        if use_pacts:
            expected = self.district_pact_expected
            results = self.district_pact_results.copy()
        else:
            expected = self.district_party_expected
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

        if use_pacts:
            return 1 / np.sum(np.square(self.pact_results / self.n_seats))

        return 1 / np.sum(np.square(self.party_results / self.n_seats))

    # quota conditions
    def quota_condition(self, district = None) -> pd.DataFrame:

        if district is None:
            expected_seats = self.party_expected * self.party_results.sum()
            results = self.party_results
        else:
            expected_seats = (
                self.district_party_expected.loc[district]
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
    def plot_quota_diff(self):

        quota_df = self.quota_condition()

        fig, ax = plt.subplots()

        ax.bar(
            self.parties,
            quota_df["seats"] - quota_df["expected"],
            color=self.colors
        )

        ax.set_xlabel("Partido")
        ax.set_ylabel("Dif. entre escaños reales y esperados")
        fig.suptitle("Disproporcionalidad entre votación y escaños")

        return fig

    def plot_parliament(self):

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
            s=100
        )
        ax.text(1, 0.1, self.n_seats, fontsize=50, ha="center")
        ax.set_axis_off()
        ax.set_xlim(0, 2.5)
        ax.legend(handles=legend_elems, title="Partidos")

        return fig

    # summary of statistics
    def summary(self):

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
        districts_meet_quota = (
            len(self.districts)
            - self.quota_condition_districts()["parties_without_quota"].sum()
        )
        
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
- Número de partidos: {len(self.parties)}
- Número efectivo de partidos: {self.effective_num_parties():.2f}
- Número de pactos: {len(self.pacts)}
- Número efectivo de pactos: {self.effective_num_parties(use_pacts=True):.2f}

Quotas:
Partidos que cumplen quota: {parties_meet_quota}/{len(self.parties)}
Distritos donde todos los partidos cumplen quota: {
    districts_meet_quota}/{len(self.districts)
}
"""
        )
        
        fig_quota = self.plot_quota_diff()
        fig_parliament = self.plot_parliament()
        plt.show()