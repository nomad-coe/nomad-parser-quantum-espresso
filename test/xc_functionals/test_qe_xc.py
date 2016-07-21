#!/usr/bin/env python
import setup_paths
from QuantumEspressoXC import translate_qe_xc_num
import sys
import logging
import re
import pprint


LOGGER = logging.getLogger(__name__)
PP = pprint.PrettyPrinter(indent=2, width=190)


def process_line(line):
    lr = line.rstrip()
    m = re.match(r".*\(([^\)]+)\)",lr)
    if m is None:
        LOGGER.error("unrecognized line: %s",lr)
        return None
    qe_xc=translate_qe_xc_num(m.group(1))
    if qe_xc is None:
        print("qe_xc: %s\nNone" % (lr))
    else:
        print(("qe_xc: %s\n" % (lr)) + PP.pprint(qe_xc))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    for line in sys.stdin:
        process_line(line)
