"""
Environmental similarity metrics (Gower, 1971; Ramirez-Villegas et al., 2011)
and deterministic Top-1 homologue matching.
"""

import math

from rufas_lactation_equivalence.biophysics import normalize_profile
from rufas_lactation_equivalence.models import (
    EquivalenceResult,
    MatchCandidate,
    NormalizedVector,
    RawEdaphoclimaticProfile,
)


def gower_similarity(vec_a: NormalizedVector, vec_b: NormalizedVector) -> float:
    """
    Computes Gower's similarity coefficient in [0.0, 1.0] across the 6 normalized dimensions.

    S_Gower(a, b) = 1 - (1/6) * sum(|a_k - b_k|)
    """
    coords_a = vec_a.as_tuple()
    coords_b = vec_b.as_tuple()
    n = len(coords_a)
    manhattan = sum(abs(a - b) for a, b in zip(coords_a, coords_b))
    return max(0.0, min(1.0, 1.0 - (manhattan / n)))


def euclidean_similarity(vec_a: NormalizedVector, vec_b: NormalizedVector) -> float:
    """
    Computes normalized Euclidean similarity in [0.0, 1.0] across the 6 normalized dimensions.

    S_Eucl(a, b) = 1 - sqrt((1/6) * sum((a_k - b_k)^2))
    """
    coords_a = vec_a.as_tuple()
    coords_b = vec_b.as_tuple()
    n = len(coords_a)
    mean_sq = sum((a - b) ** 2 for a, b in zip(coords_a, coords_b)) / n
    return max(0.0, min(1.0, 1.0 - math.sqrt(mean_sq)))


def compute_all_equivalences(
    brazil_profiles: list[RawEdaphoclimaticProfile],
    us_profiles: list[RawEdaphoclimaticProfile],
    region_adjustments: dict[str, dict[str, float]],
    top_k: int = 3,
) -> list[EquivalenceResult]:
    """
    Computes the normalized vectors and ranks candidate US states for every Brazilian UF.
    Only US states mapped to an active RuFaS region (not None and != "None") are eligible candidates.
    """
    valid_us: list[tuple[RawEdaphoclimaticProfile, NormalizedVector]] = []
    for us_prof in us_profiles:
        if us_prof.rufas_region and us_prof.rufas_region != "None":
            if us_prof.rufas_region not in region_adjustments:
                raise ValueError(
                    f"US state {us_prof.code} references unknown RuFaS region '{us_prof.rufas_region}'"
                )
            valid_us.append((us_prof, normalize_profile(us_prof)))

    if not valid_us:
        raise ValueError("No valid US candidate states provided.")

    results: list[EquivalenceResult] = []
    for br_prof in sorted(brazil_profiles, key=lambda p: int(p.code)):
        br_vec = normalize_profile(br_prof)
        candidates: list[MatchCandidate] = []
        for us_prof, us_vec in valid_us:
            g_sim = gower_similarity(br_vec, us_vec)
            e_sim = euclidean_similarity(br_vec, us_vec)
            candidates.append(
                MatchCandidate(
                    us_code=us_prof.code,
                    us_abbrev=us_prof.abbrev,
                    us_name=us_prof.name,
                    rufas_region=str(us_prof.rufas_region),
                    gower_similarity=g_sim,
                    euclidean_similarity=e_sim,
                )
            )

        candidates.sort(
            key=lambda c: (
                round(c.gower_similarity, 12),
                round(c.euclidean_similarity, 12),
                -int(c.us_code),
            ),
            reverse=True,
        )

        best = candidates[0]
        adjustments = dict(region_adjustments[best.rufas_region])
        results.append(
            EquivalenceResult(
                brazil_profile=br_prof,
                brazil_vector=br_vec,
                best_match=best,
                top_candidates=tuple(candidates[: max(1, top_k)]),
                wood_adjustments=adjustments,
            )
        )

    return results
