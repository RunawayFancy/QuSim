import numpy as np
from typing import Dict, List, Tuple, Optional

from qusim.pulse_gen.noise_gen import noise_gen, clear_corr_noise_cache

Channel = Tuple[str, object]  # e.g. ("INT", [0,1]) or ("XY", 0)


def build_channel_noise_series(
    channel_noise: List[Tuple[Channel, list]],
    simopt,
    waveform_map: Optional[Dict[Channel, np.ndarray]] = None,
    reset_cache: bool = True,
) -> Dict[Channel, Dict[str, np.ndarray]]:
    """
    Build per-channel noise arrays using the same path as runtime noise_gen.

    Returns:
        {channel: {"raw": raw_noise, "applied": applied_noise}}
    where raw_noise is the sum of per-config noise (before sum/multiply),
    and applied_noise is the effective additive term applied to a waveform.
    """
    if reset_cache:
        clear_corr_noise_cache()

    tlist = simopt.tlist
    zero_wf = np.zeros_like(tlist, dtype="float64")
    results: Dict[Channel, Dict[str, np.ndarray]] = {}

    for chan, noise_cfg_list in channel_noise:
        wf = zero_wf if waveform_map is None else waveform_map.get(chan, zero_wf)
        raw_sum = np.zeros_like(tlist, dtype="float64")
        applied = np.zeros_like(tlist, dtype="float64")

        for cfg in noise_cfg_list:
            n = np.real(noise_gen(cfg, wf.copy()))
            raw_sum += n
            if cfg.methods == "sum":
                applied += n
            elif cfg.methods == "multiply":
                applied *= n
            else:
                raise AttributeError("noise config missing methods or invalid methods")
        # print(chan)
        # print(results)
        results[str(chan)] = {"raw": raw_sum, "applied": applied}

    return results


def correlation_report(
    channel_noise: List[Tuple[Channel, list]],
    simopt,
    corr_id: str,
    waveform_map: Optional[Dict[Channel, np.ndarray]] = None,
    reset_cache: bool = True,
) -> Dict[str, object]:
    """
    Compute pairwise Pearson correlations among channels that share corr_id.
    """
    series = build_channel_noise_series(
        channel_noise=channel_noise,
        simopt=simopt,
        waveform_map=waveform_map,
        reset_cache=reset_cache,
    )

    channels = []
    for chan, cfgs in channel_noise:
        if any(getattr(cfg, "corr_id", None) == corr_id for cfg in cfgs):
            channels.append(chan)

    corr_mat = np.eye(len(channels))
    for i in range(len(channels)):
        for j in range(i + 1, len(channels)):
            a = series[str(channels[i])]["raw"]
            b = series[str(channels[j])]["raw"]
            if np.std(a) == 0 or np.std(b) == 0:
                c = 0.0
            else:
                c = float(np.corrcoef(a, b)[0, 1])
            corr_mat[i, j] = c
            corr_mat[j, i] = c

    return {
        "corr_id": corr_id,
        "channels": channels,
        "corr_matrix": corr_mat,
    }


def assert_correlated(
    channel_noise: List[Tuple[Channel, list]],
    simopt,
    corr_id: str,
    min_corr: float = 0.99,
    waveform_map: Optional[Dict[Channel, np.ndarray]] = None,
    reset_cache: bool = True,
) -> bool:
    """
    Return True if all channel pairs with corr_id have correlation >= min_corr.
    """
    report = correlation_report(
        channel_noise=channel_noise,
        simopt=simopt,
        corr_id=corr_id,
        waveform_map=waveform_map,
        reset_cache=reset_cache,
    )
    cm = report["corr_matrix"]
    if cm.size == 1:
        return True
    tri = cm[np.triu_indices_from(cm, k=1)]
    return bool(np.all(tri >= min_corr))
