import setup_paths
import re
import sys
import os
import logging
from nomadcore.match_highlighter import ANSI

LOGGER = logging.getLogger(__name__)


# regex for _valid_ fortran float output, what a mess ...
RE_f = (r"(?:" + # outer alternative, between numbers and number-too-wide-for-field markers
    r"([+-]?)(?:" + # MANTISSA, SIGN (group 1, optional), followed by alternatives
        '|'.join([ # MANTISSA
            r"(\d+(?!\.))", # MANTISSA without a decimal point, group 2
            r"(\d*)" + ( # MANTISSA, WHOLE part (group 3)
                # we need negative look-ahead/look-behind assertions around the
                # decimal point as there is too much optional stuff around
                r"(?<![^\d\s+-])" + # char preceding the dot must be nothing but number, whitespace, or sign
                r"\." +
                r"(?![^eEdD\d\s,])" + # char succeeding the dot must be nothing but number, exponential/precision char, comma or whitespace
                r"(\d*)" # MANTISSA, FRACTIONAL part (group 4), separated by dot
            )
        ]) +
    r")(?:" + ( # EXPONENT part (optional)
        r"([eEdD])" + # PRECISION (group5)
        r"([+-]?)(\d*)" # EXPONENT SIGN (group 6), VALUE (group 7)
    ) + ")?" + # make precision/exponet part optinal
    r"|(\*+))" # outer alternative, between numbers and number-too-wide markers (group 8)
)
cRE_f = re.compile(RE_f)


def match_to_float(m, group_offset=0):
    group = [ m.group(0) ] + [ m.group(group_offset + i) for i in range(1,9)]
    LOGGER.debug("g: %s", str(group))
    if group[8] is not None:
        pyfloat_str = 'nan'
        dtype = 'f'
    else:
        pyfloat_str = group[1] # sign, maybe zero-length
        if group[2] is not None:
            pyfloat_str += group[2]
            dtype = 'i'
        else:
            pyfloat_str += group[3] if len(group[3])>0 else '0'
            pyfloat_str += '.'
            pyfloat_str += group[4] if len(group[4])>0 else '0'
            dtype = 'f'
        if group[5] is not None:
            pyfloat_str += 'e' + group[6]
            pyfloat_str += group[7] if len(group[7])>0 else '0'
            dtype = 'f'
    LOGGER.debug("pyfloat_str: %s", pyfloat_str)
    return (float(pyfloat_str), dtype)


RE_unescape = {
    '"': re.compile(r'""'),
    "'": re.compile(r"''"),
}


def unquote_string(value):
    result = value[1:-1]
    result = RE_unescape[value[0]].sub(value[0], result)
    return result


# quoted strings
cRE_string_quoted = re.compile(r"(?:'[^']*'|\"[^\"]*\")")
cRE_comment = re.compile(r"\s*!.*")
RE_identifier = r"[a-zA-Z]\w*" # fortran identifier
cRE_start_group = re.compile(r'\s*&(' + RE_identifier + r')') # beginning of namelist group 
cRE_end_group = re.compile(r'\s*/')
cRE_start_assignment = re.compile(r'\s*(?P<target>' + RE_identifier + r')(?:\(\s*(?P<subscript>[^\)]*?)\s*\))?\s*=\s*')
cRE_assigned_value = re.compile(
    r'\s*(?:' + '|'.join([
        r'(?P<num>' + RE_f + r')', # integers and floats
        r'\(\s*(?P<cnum_r>' + RE_f + r')\s*,\s*(?P<cnum_i>' + RE_f + r')\s*\)', # complex numbers
        r'(?P<bool_t>\.t(?:rue)?\.)', # true-value bool
        r'(?P<bool_f>\.f(?:alse)?\.)', # false-value bool
        r"(?P<str_s>'[^']*(?:[^']|'')*'(?!'))", # single-quoted string, closed, allowing for escaped quotes ('')
        r'(?P<str_d>"[^"]*(?:[^"]|"")*"(?!"))', # double-quoted string, closed, allowing for escaped quotes ("")
        r"(?P<str_s_nc>'[^']*(?:[^']|'')*)", # single-quoted string, not closed
        r'(?P<str_d_nc>"[^"]*(?:[^"]|"")*)', # double-quoted string, not closed
        r'(?P<comment>!.*)', # comment
    ]) + ')', re.I)
cRE_str_s_close = re.compile(r"([^']*(?:[^']|'')*'(?!'))") # single-quoted string, closing
cRE_str_d_close = re.compile(r'([^"]*(?:[^"]|"")*"(?!"))') # double-quoted string, closing
cRE_comma = re.compile(r'\s*,')


class FortranNamelistParser(object):
    """Parser for Fortran 90 Namelists
    """
    def __init__(self, file_path):
        self.input_tree = {}
        self.file_path = file_path
        self.state = 0
        self.__nl_group = None
        self.__target = None
        self.__subscript = None
        self.__values = None
        self.__types = None
        self.__nvalues_after_comma = 0

    def parse(self):
        with open(self.file_path, "r") as fIn:
            # split lines into 'line' and 'comment' parts
            for line in fIn:
                # strip final newline if it exists
                if line[-1] == '\n':
                    line = line[:-1]
                self.parse_line(line)

    def parse_line(self, line):
        last_end = 0
        while last_end<len(line):
            if self.state == 0:
                # we have no open group
                m = cRE_start_group.match(line, last_end)
                if m is not None:
                    self.__nl_group = m.group(1)
                    sys.stdout.write(ANSI.FG_BRIGHT_GREEN + m.group() + ANSI.RESET)
                    last_end = m.end()
                    self.state = 1
                    self.onOpen_namelist_group(m.group(1))
                    continue
                # but comments may appear here
                m = cRE_comment.match(line, last_end)
                if m is not None:
                    sys.stdout.write(ANSI.FG_BLUE + m.group() + ANSI.RESET)
                    last_end = m.end()
                    self.onComment(m.group())
                    continue
            elif self.state==3:
                # we are inside single-quoted multiline string
                m = cRE_str_s_close.match(line, last_end)
                if m is None:
                    sys.stdout.write(ANSI.FG_YELLOW + line[last_end:] + ANSI.RESET)
                    self.__values[-1] += "\n" + line
                    last_end=len(line)
                else:
                    sys.stdout.write(ANSI.FG_YELLOW + m.group() + ANSI.RESET)
                    self.__values[-1] += "\n" + m.group(1)
                    self.__values[-1] = unquote_string(self.__values[-1])
                    self.__types[-1] = 'C'
                    last_end=m.end()
                    self.state = 2
                    continue
            elif self.state==4:
                # we are inside double-quoted multiline string
                m = cRE_str_d_close.match(line, last_end)
                if m is None:
                    sys.stdout.write(ANSI.FG_YELLOW + line[last_end:] + ANSI.RESET)
                    self.__values[-1] += "\n" + line
                    last_end=len(line)
                else:
                    sys.stdout.write(ANSI.FG_YELLOW + m.group() + ANSI.RESET)
                    self.__values[-1] += "\n" + m.group(1)
                    self.__values[-1] = unquote_string(self.__values[-1])
                    self.__types[-1] = 'C'
                    last_end=m.end()
                    self.state = 2
                    continue
            elif self.state == 1 or self.state == 2:
                # we are inside opened group
                #   check for group-closing /
                m = cRE_end_group.match(line, last_end)
                if m is not None:
                    if self.__target is not None:
                        self.onClose_value_assignment(
                            self.__target, self.__subscript,
                            self.__values, self.__types)
                    self.__target = None
                    self.__subscript = None
                    self.__values = None
                    self.__types = None
                    self.__nvalues_after_comma = 0
                    self.onClose_namelist_group(self.__nl_group)
                    self.__nl_group = None
                    sys.stdout.write(ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_GREEN + m.group() + ANSI.RESET)
                    self.state = 0
                    last_end = m.end()
                    continue
                #   check for new assignment
                m = cRE_start_assignment.match(line, last_end)
                if m is not None:
                    if self.__target is not None:
                        self.onClose_value_assignment(
                            self.__target, self.__subscript,
                            self.__values, self.__types)
                    self.state = 2
                    last_end=m.end()
                    sys.stdout.write(ANSI.FG_GREEN + m.group() + ANSI.RESET)
                    self.__target = m.group('target')
                    self.__subscript = m.group('subscript')
                    self.__values = []
                    self.__types = []
                    self.__nvalues_after_comma = 0
                    self.onOpen_value_assignment(
                        self.__target, self.__subscript)
                    continue
                if self.state == 2:
                    # we are inside the values-part of an assignment
                    m = cRE_assigned_value.match(line, last_end)
                    if m is not None:
                        if m.group('num') is not None:
                            (value, dtype) = match_to_float(m, group_offset=1)
                            self.__values.append(value)
                            self.__types.append(dtype)
                        elif m.group('cnum_r') is not None:
                            (cnum_r, dtype) = match_to_float(m, group_offset=10)
                            (cnum_i, dtype) = match_to_float(m, group_offset=19)
                            self.__values.append(complex(cnum_r, cnum_i))
                            self.__types.append('complex')
                        elif m.group('bool_t') is not None:
                            self.__values.append(True)
                            self.__types.append('b')
                        elif m.group('bool_f') is not None:
                            self.__values.append(False)
                            self.__types.append('b')
                        elif m.group('str_s') is not None:
                            self.__values.append(unquote_string(m.group('str_s')))
                            self.__types.append('C')
                        elif m.group('str_d') is not None:
                            self.__values.append(unquote_string(m.group('str_d')))
                            self.__types.append('C')
                        elif m.group('str_s_nc') is not None:
                            # non-closed single-quoted string
                            self.state=3
                            self.__values.append(m.group('str_s_nc'))
                            self.__types.append('string_singlequoted')
                        elif m.group('str_d_nc') is not None:
                            # non-closed double-quoted string
                            self.state=4
                            self.__values.append(m.group('str_d_nc'))
                            self.__types.append('string_doublequoted')
                        elif m.group('comment') is not None:
                            sys.stdout.write(ANSI.FG_BLUE + m.group() + ANSI.RESET)
                            last_end=m.end()
                            self.onComment(m.group())
                            continue
                        self.__nvalues_after_comma +=1
                        sys.stdout.write(ANSI.FG_YELLOW + m.group() + ANSI.RESET)
                        last_end=m.end()
                        continue
                    # special meaning of comma: may indicate Null values in array assignments
                    m = cRE_comma.match(line, last_end)
                    if m is not None:
                        if self.__nvalues_after_comma is 0:
                            self.__values.append(None)
                            self.__types.append(None)
                        self.__nvalues_after_comma = 0
                        sys.stdout.write(ANSI.FG_MAGENTA + m.group() + ANSI.RESET)
                        last_end = m.end()
                        continue
            break
        if last_end < len(line):
            if self.state > 0 and self.state < 5 and line[last_end:].strip():
                # states we as the base class are handling, but with leftover chars on a line
                LOGGER.error("ERROR: leftover chars in line while inside namelist group: '%s'", line[last_end:])
                sys.stdout.write(ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_RED + line[last_end:] + ANSI.RESET)
            else:
                sys.stdout.write(ANSI.BEGIN_INVERT + ANSI.FG_BLUE + line[last_end:] + ANSI.RESET)
        sys.stdout.write('\n')

    # Hooks to be overloaded in derived classes in order to do stuff
    def onComment(self, comment):
        pass

    def onOpen_namelist_group(self, groupname):
        pass

    def onClose_namelist_group(self, groupname):
        pass

    def onOpen_value_assignment(self, target, subscript):
        pass

    def onClose_value_assignment(self, target, subscript, values, dtypes):
        if subscript is None:
            LOGGER.error("SET %s = %s (types: %s)", target, str(values), str(dtypes))
        else:
            LOGGER.error("SET %s(%s) = %s (types: %s)", target, subscript, str(values), str(dtypes))

if __name__ == "__main__":
    parser = FortranNamelistParser(sys.argv[1])
    parser.parse()
