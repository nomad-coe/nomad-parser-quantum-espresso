#!/usr/bin/env python
import setup_paths
from QuantumEspressoXC import translate_qe_xc_num
import sys
import logging
import re
import json

LOGGER = logging.getLogger(__name__)

re_line = re.compile(
    r"\s*\d+"
    r"\s+(?P<exact_exchange_fraction>[0-9\.]+)"
    r"\s+.*\((?P<xc_functional_num>[^\)]+)\)"
)

re_comment = re.compile(r"^\s*#.*")

def process_line(line):
    if re_comment.match(line):
        return
    lr = line.rstrip()
    print("qe_xc: %s" % (lr))
    m = re_line.match(lr)
    if m is None:
        LOGGER.error("unrecognized line:",lr)
        return None
    qe_xc = None
    try:
        qe_xc=translate_qe_xc_num(
            m.group('xc_functional_num'),
            float(m.group('exact_exchange_fraction'))
        )
    except RuntimeError as e:
        print("  Error: %s" % (str(e)))
    if qe_xc is None:
        print("  None")
    else:
        print(json.dumps(qe_xc, indent=2, sort_keys=True))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    for line in sys.stdin:
        process_line(line)
