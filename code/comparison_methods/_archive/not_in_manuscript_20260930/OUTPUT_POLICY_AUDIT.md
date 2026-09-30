# Native-output constraint: targeted audit, 2026-09-28

The user permits direct native API calls and our timing, but prohibits extending
or changing model outputs. Missing capabilities must be reported as unavailable.
This is a targeted audit, not certification of every historical comparison file.

## Compute

The replacement scripts in `compute/*_native.py` contain only direct public API
calls, ordinary input preparation, read-only checks, and performance reporting.
No spot compositions are computed. scCube coordinates and SPIDER's expression
view remain exactly as returned. Earlier custom-wrapper measurements are excluded.

## Existing overview export

`overview/generate.py` previously reconstructed SPIDER circular-capture
composition as membership multiplied by one-hot cell labels. The installed
st-spider 1.2.0 circular API does not return that composition. The fallback has
been removed. The legacy composition-dependent export now stops with an explicit
unavailable-output error for that branch. It must not silently substitute a
computed composition or label it native. Already saved derived fields are not
made native by this code change and must not be presented as native capabilities.

The square branch returns its count matrix as the fourth native return value.
scCube separately provides an upstream `calculate_spot_prop` API; invoking that
upstream API is different from reconstructing a missing SPIDER output ourselves.
Any such optional native call must be named in provenance. Display summaries
such as a dominant-type color must remain separate from unchanged native outputs
and must not be described as native labels or new ground truth.

Historical adapter-based exports also normalize schemas, coordinates, and labels.
Those have not been regenerated or certified by the compute repair. Do not reuse
them as proof that a method natively provides the standardized fields.
