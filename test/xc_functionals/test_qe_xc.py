#!/usr/bin/env python3
#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import setup_paths
from QuantumEspressoXC import translate_qe_xc_num
import sys
import logging
import re
import json

LOGGER = logging.getLogger(__name__)

re_line = re.compile(
    r"^\s*\d+\s+(?P<testcase>" + (
        r"(?P<exact_exchange_fraction>[0-9\.]+)"
        r"\s+.*\((?P<xc_functional_num>[^\)]+)\)"
    ) + r")\s*$"
)

re_comment = re.compile(r"^\s*#.*")

def process_line(line):
    if re_comment.match(line):
        return
    m = re_line.match(line)
    if m is None:
        LOGGER.error("unrecognized line:",line)
        return None
    print("qe_xc: %s" % (m.group('testcase')))
    qe_xc = None
    try:
        (method_qe_xc, qe_xc) = translate_qe_xc_num(
            m.group('xc_functional_num'),
            float(m.group('exact_exchange_fraction'))
        )
    except RuntimeError as e:
        print("  Error: %s" % (str(e)))
    if qe_xc is None:
        print("  None")
    else:
        print(json.dumps(method_qe_xc, indent=2, sort_keys=True))
        print(json.dumps(qe_xc, indent=2, sort_keys=True))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    for line in sys.stdin:
        process_line(line)
