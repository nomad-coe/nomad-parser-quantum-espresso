#!/usr/bin/env python
import setup_paths
from QuantumEspressoXC import translate_qe_xc_num
import sys
import logging
import re


LOGGER = logging.getLogger(__name__)


def process_line(line):
    lr = line.rstrip()
    m = re.match(r".*\(([^\)]+)\)",lr)
    if m:
        translate_qe_xc_num(m.group(1))
    else:
        LOGGER.error("unrecognized line: %s",lr)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    for line in sys.stdin:
        process_line(line)
