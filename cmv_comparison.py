"""Compare independent CMV codon evidence without promoting review observations.

PASS here means only that configured exploratory read-level thresholds were met.
"""

PASS_STATUS = 'READ_PHASED_EXPLORATORY'


def compare_calls(legacy_calls, read_calls, panel_mutations=()):
    legacy = {(v.get('gene'), v.get('mutation')) for v in legacy_calls if v.get('gene') and v.get('mutation')}
    panel = set(panel_mutations)
    # One biological substitution may be represented by multiple codon paths.
    by_key = {}
    for v in read_calls:
        key = (v.get('gene'), v.get('mutation'))
        if not all(key):
            continue
        by_key.setdefault(key, []).append(v)

    summary = []
    for key in sorted(legacy | panel | {k for k, vs in by_key.items() if any(v.get('status') == PASS_STATUS for v in vs)}):
        observations = by_key.get(key, [])
        passed = [v for v in observations if v.get('status') == PASS_STATUS]
        best = max(passed or observations, key=lambda v: (int(v.get('alt_reads') or 0), int(v.get('depth') or 0)), default={})
        if key in legacy and passed:
            comparison = 'BOTH_CALLERS'
        elif passed:
            comparison = 'READ_LEVEL_ONLY_PASS_REVIEW'
        elif key in legacy and observations:
            comparison = 'LEGACY_WITH_READ_LEVEL_REVIEW'
        elif key in legacy:
            comparison = 'LEGACY_ONLY'
        elif observations:
            comparison = 'PANEL_READ_LEVEL_REVIEW'
        else:
            comparison = 'PANEL_NOT_OBSERVED'
        summary.append({**best, 'gene': key[0], 'mutation': key[1], 'comparison': comparison,
                        'read_level_pass': bool(passed), 'legacy_present': key in legacy,
                        'read_level_observations': len(observations)})

    # Full evidence remains separate and never acquires a misleading comparison status.
    details = []
    for key in sorted(by_key):
        for v in sorted(by_key[key], key=lambda x: (x.get('alternate_codon') or '', -(x.get('alt_reads') or 0))):
            details.append({**v, 'legacy_present': key in legacy,
                            'read_level_pass': v.get('status') == PASS_STATUS})
    return summary, details
