"""Read SEAKMC discovery records; NEB separately tests prescribed paths."""
import csv
import math
from pathlib import Path


def seakmc_metrics(directory):
    directory = Path(directory)
    summary = directory / 'Seakmc_summary.csv'
    paths = [summary] if summary.is_file() else sorted((directory / 'SPOut').glob('KMC_*_SPs.csv'))
    forward, backward, counts = [], [], []
    records = 0
    for path in paths:
        with path.open(encoding='utf-8', errors='replace') as handle:
            for row in csv.DictReader(handle):
                records += 1
                normalized = {key.lower().replace('_', ''): value for key, value in row.items() if key}
                for target, names in ((forward, ('forwardbarrier', 'barrier')), (backward, ('backwardbarrier', 'backbarrier')), (counts, ('nsps',))):
                    for name in names:
                        try:
                            value = float(normalized[name])
                        except (KeyError, ValueError, TypeError):
                            continue
                        if math.isfinite(value) and value >= 0:
                            target.append(value)
                            break
    result = {}
    if forward:
        result['min_forward_barrier'] = min(forward)
    if backward:
        result['min_backward_barrier'] = min(backward)
    if records:
        result['saddle_records'] = records
    if counts:
        result['reported_search_saddle_count'] = sum(counts)
    return result
