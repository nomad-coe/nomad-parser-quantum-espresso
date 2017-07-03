import setup_paths
import calendar
import json
import os
import re
import numpy as np
import logging
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM, CachingLevel

############################################################
# This file contains functions and constants that are needed
# by more than one parser.
############################################################


LOGGER = logging.getLogger(__name__)

# fortran float, alternate too-long-for-field fortran marker
RE_f = r"(?:[+-]?\d+(?:\.\d+)?(?:[eEdD][+-]?\d+)?|\*+)"
cRE_f = re.compile(RE_f)
# fortran int, alternate too-long-for-field fortran marker
RE_i = r"(?:[+-]?\d+|\*+)"
cRE_i = re.compile(RE_i)
NAN = float('nan')

def re_vec(name, units='', split="\s+"):
    """generator for 3-component vector regex"""
    if units:
        units = '__' + units
    res = (
        r'(?P<' + name + r'_x' + units + r'>' + RE_f + r')' + split +
        r'(?P<' + name + r'_y' + units + r'>' + RE_f + r')' + split +
        r'(?P<' + name + r'_z' + units + r'>' + RE_f + r')'
        )
    return res


# loading metadata from
# nomad-meta-info/meta_info/nomad_meta_info/quantum_espresso.nomadmetainfo.json
META_INFO = loadJsonFile(
    filePath=os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../../../nomad-meta-info/meta_info/nomad_meta_info/quantum_espresso.nomadmetainfo.json")),
    dependencyLoader=None,
    extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
    uri=None)[0]

PARSER_INFO_DEFAULT = {
  "name": "parser_quantum_espresso",
  "version": "0.0.1"
}

# constants for date conversion
MONTHS = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
           'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ]
MONTH_NUMBER = { MONTHS[num]: num+1 for num in range(0,12) }


QE_SMEARING_KIND = {
    '-99': 'fermi',
     '-1': 'marzari-vanderbilt',
      '0': 'gaussian',
      '1': 'methfessel-paxton',
    'Marzari-Vanderbilt smearing': 'marzari-vanderbilt',
     'Methfessel-Paxton smearing': 'methfessel-paxton',
              'gaussian smearing': 'gaussian',
           'Fermi-Dirac smearing': 'fermi',
             'tetrahedron method': 'tetrahedra',
}


class ParserQuantumEspresso(object):
    """Base class for all Quantum Espresso parsers"""
    def __init__(self,cachingLevelForMetaName=None, coverageIgnoreList=None,
                 re_program_name=None):
        self.re_program_name = re_program_name
        self.parserInfo = PARSER_INFO_DEFAULT.copy()
        self.cachingLevelForMetaName = {}
        for name in META_INFO.infoKinds:
            # set all temporaries to caching-only
            if name.startswith('x_qe_t_'):
                self.cachingLevelForMetaName[name] = CachingLevel.Cache
        # common prosa in espresso output
        self.coverageIgnoreList = [
            # ignore empty lines
            r"\s*",
            # table separators
            r"^\s*[=%-]+\s*$",
            r"^\s*%\s*%\s*$",
        ]
        self.coverageIgnore = None

    def parse(self):
        self.coverageIgnore = re.compile(r"^(?:" + r"|".join(self.coverageIgnoreList) + r")$")
        mainFunction(self.mainFileDescription(), META_INFO, self.parserInfo,
                    cachingLevelForMetaName=self.cachingLevelForMetaName,
                    superContext=self)

    def adHoc_suicide_qe_program_name(self, parser):
        if self.re_program_name is not None:
            if not self.re_program_name.match(
                    parser.lastMatch['x_qe_program_name']):
                raise Exception(
                    "mainFile program name was: %s, unsuited for %s" % (
                        parser.lastMatch['x_qe_program_name'],
                        type(self).__name__))

    def mainFileDescription(self):
        # assemble matchers and submatchers
        result = SM(
            name='root',
            weak=True,
            startReStr="",
            subMatchers=[
                SM(name='ktab_cIgn', coverageIgnore=True,
                   # early output seen in benchmark.out.v5.3.0.inp\=vdw1.in.1452257026
                   startReStr=r"\s*Generating kernel table - May take several minutes.*$",
                ),
                SM(name='newRun', repeats=True, required=True,
                   startReStr=(
                       # program name, e.g. PWSCF, DOS
                       r"\s*Program\s+(?P<x_qe_program_name>\S+)\s+v\." +
                       # version
                       r"(?P<program_version>\S+(?:\s+\(svn\s+rev\.\s+\d+\s*\))?)" +
                       r"\s+starts" +
                       # newer espresso: "on $date"
                       r"(?:(?:\s+on\s+(?P<time_run_date_start__strQeDate>.+?)?)\s*$|" +
                       # older espresso has just "..." and date on new line
                       r"(?:\s*\.\.\.)\s*$)"
                   ),
                   adHoc = self.adHoc_suicide_qe_program_name,
                   fixedStartValues={'program_name': 'Quantum Espresso',
                                     'program_basis_set_type': 'plane waves'},
                   sections=['section_run'],
                   subMatchers=([
                       # older espresso versions have start date on separate line
                       SM(name='run_date',
                          startReStr=r"\s*Today is\s*(?P<time_run_date_start__strQeDate>.+?)\s*$"
                       ),
                       SM(name='copyright_msg', coverageIgnore=True,
                          # ignore copyright/citation msg
                          startReStr=r"\s*This program is part of the open-source Quantum ESPRESSO suite",
                          subMatchers=[
                              SM(name='copyright_msg010', coverageIgnore=True,
                                 startReStr=r"\s*for quantum simulation of materials; please (?:cite|acknowledge)",
                              ),
                              SM(name='copyright_msg020', coverageIgnore=True,
                                 startReStr=r"\s*\"P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 \(2009\);\s*",
                              ),
                              SM(name='copyright_msg030', coverageIgnore=True,
                                 startReStr=r"\s*URL http://www.quantum-espresso.org\",\s*",
                              ),
                              SM(name='copyright_msg040', coverageIgnore=True,
                                 startReStr=r"\s*in publications or presentations arising from this work. More details at",
                              ),
                              SM(name='copyright_msg050', coverageIgnore=True,
                                 startReStr=r"\s*http://www.quantum-espresso.org/quote(?:\.php)?",
                              ),
                              SM(name='copyright_msg055', coverageIgnore=True,
                                 startReStr=r"\s*http://www.quantum-espresso.org/wiki/index.php/Citing_Quantum-ESPRESSO\s*",
                              ),
                          ],
                       ),
                   ] + self.run_submatchers() + [
                       SM(name='end_date',
                          startReStr=r"\s*This run was terminated on:\s*(?P<time_run_date_end__strQeDate>.+?)\s*$",
                       ),
                       SM(name='job_done',
                          startReStr=r"\s*JOB DONE\.\s*",
                          adHoc=lambda p: p.backend.addValue('run_clean_end', True),
                       ),
                   ]),
                )
            ]
        )
        return result

    def run_submatchers(self):
        return []

    def strValueTransform_strQeDate(self, espresso_date):
        if espresso_date is None:
            return None
        epoch = 0
        match = re.match(
            r"(\d+)\s*([A-Za-z]+)\s*(\d+)\s+at\s+(\d+):\s*(\d+):\s*(\d+)",
            espresso_date)
        if match:
            month = MONTH_NUMBER[match.group(2)]
            epoch = calendar.timegm(
                (int(match.group(3)), int(month), int(match.group(1)),
                 int(match.group(4)), int(match.group(5)), int(match.group(6))))
        else:
            match = re.match(
                r"\s*(\d+):\s*(\d+):\s*(\d+)\s+(\d+)\s*([A-Za-z]+)\s*(\d+)\s*",
                espresso_date)
            if match:
                month = MONTH_NUMBER[match.group(5)]
                epoch = calendar.timegm(
                    (int(match.group(6)), int(month), int(match.group(4)),
                     int(match.group(1)), int(match.group(2)), int(match.group(3))))
            else:
                raise RuntimeError("unparsable date: %s", espresso_date)
        return(epoch)
    strValueTransform_strQeDate.units = 's'

    def strValueTransform_strQeTimespan(self, espresso_timespan):
        if espresso_timespan is None:
            return None
        match = re.match(r"(?:\s*(\d+)\s*m)?" + # minutes (optional)
                         r"(?:\s*(" + RE_f + ")\s*s)?", # seconds (optional)
                         espresso_timespan)
        if not match:
            raise RuntimeError(
                "unparsable timespan (regex match failed): %s",
                espresso_timespan)
        had_time_components = 0
        timespan_seconds = 0.0
        if match.group(2) is not None:
            timespan_seconds += float(match.group(2))
            had_time_components += 1
        if match.group(1) is not None:
            timespan_seconds += float(match.group(1))*60
            had_time_components += 1
        if had_time_components == 0:
            raise RuntimeError(
                "unparsable timespan (no time components extracted): %s",
                espresso_timespan)
        return(timespan_seconds)
    strValueTransform_strQeTimespan.units = 's'

    def addSectionDict(self, backend, section_name, section_dict):
        gIndex = backend.openSection(section_name)
        for key, value in sorted(section_dict.items()):
            if isinstance(value, (list,dict)):
                backend.addValue(key, value)
            else:
                backend.addValue(key, value)
        backend.closeSection(section_name, gIndex)
