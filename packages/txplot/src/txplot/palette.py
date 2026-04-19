"""Color palettes. Defaults are colorblind-safe (Okabe-Ito); larger sets
delegate to glasbey."""
from __future__ import annotations

__all__ = ["OKABE_ITO", "CONSEQUENCE_COLORS", "categorical"]

OKABE_ITO = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # gray
]

CONSEQUENCE_COLORS: dict[str, str] = {
    "missense": "#D55E00",
    "nonsense": "#000000",
    "synonymous": "#56B4E9",
    "splice": "#CC79A7",
    "frameshift": "#D55E00",
    "inframe_indel": "#E69F00",
    "start_lost": "#000000",
    "stop_lost": "#000000",
    "stop_gained": "#000000",
    "5_prime_utr": "#999999",
    "3_prime_utr": "#999999",
    "intronic": "#BBBBBB",
    "other": "#666666",
}


def categorical(n: int, colorblind_safe: bool = True) -> list[str]:
    """Return a list of ``n`` distinct hex colors.

    For ``n <= 8``, uses Okabe-Ito. For larger sets, delegates to glasbey.
    """
    if n <= len(OKABE_ITO):
        return OKABE_ITO[:n]
    import glasbey  # lazy; only needed for large palettes

    return glasbey.create_palette(palette_size=n, colorblind_safe=colorblind_safe, as_hex=True)
