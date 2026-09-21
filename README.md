# Named AgNP observation-channel scores as forcing provenance under evidence gates

**Thesis #34.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #8 (papaya AgNP observation channel) and Thesis #19 (forcing admission under evidence gates).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

When a named papaya AgNP observation channel produces scores that are barred from Θ, which of those scores—if any—may be admitted only as provenance of a protocol-constant tip-ODE forcing without gate bypass or soft-prior leakage?

One record may. It is the class-B summary of `C_amyl` from Thesis #8. Twenty-two ledger calls return refused. One call, A01, admits that summary as provenance of `u_tip = 1`. The rank, the leading eigenvalue, the condition number, and `K_c` are not copied. Companion channels are refused because they are not the named AgNP channel. The SHA-256 of kinetic Θ is `1ace2bb9eed19f85046f25ab690076d0d5920cf19219ddb362c775be9f41d0e8` before the refusals and the same string after the admission. That string is the Thesis #19 kinetic digest. The forcing digest changes. A soft-prior contrast that would have set `d0` to 0.195683 is computed beside the ledger and is not written.

No plate table from the 2022 wet-lab thesis is re-tabulated. No Fisher matrix is recomputed. An admitted schedule is not a dose.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/admission.py` | AgNP admission ledger and inherited tip-demonstration field |
| `sim/candidates.json` | Transcribed Thesis #8 channel records (not recomputed) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Eligibility, trajectories, digests, equilibria |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/admission.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the ledger. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`. The candidate-file SHA-256 is pinned inside `sim/admission.py`. A byte edit that does not update the pin raises.

## Cite

Ogbonna KE. Named AgNP observation-channel scores as forcing provenance under evidence gates [Internet]. Thesis #34 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-34-agnp-channel-forcing-provenance

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
