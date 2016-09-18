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
cRE_start_assignment = re.compile(r'\s*(?P<target>' + RE_identifier + r')(?:\((?P<subscript>[^\)]*)\))?\s*=\s*')
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

cRE_subscript = re.compile(r'\s*,?\s*(?:(\d*):\s*(\d*)|(\d*))')


class FortranNamelistParser(object):
    """Parser for Fortran 90 Namelists
    """
    def __init__(self, file_path):
        self.input_tree = {}
        self.file_path = file_path
        self.state = 0
        self.annotateFile = sys.stdout
        self.__nl_group = None
        self.__target = None
        self.__subscript = None
        self.__values = None
        self.__types = None
        self.__nvalues_after_comma = 0
        self.bad_input = False

    def parse(self):
        with open(self.file_path, "r") as fIn:
            # split lines into 'line' and 'comment' parts
            for line in fIn:
                # strip final newline if it exists
                #   (may not exist at end of file)
                if line[-1] == '\n':
                    line = line[:-1]
                self.parse_line(line)
        if self.bad_input:
            self.onBad_input()

    def parse_subscript(self, subscript):
        if subscript is None:
            return None
        result = []
        last_end = 0
        while last_end<len(subscript):
            m = cRE_subscript.match(subscript, last_end)
            if m is None:
                break
            elif m.group(3) is not None:
                if self.annotateFile:
                    self.annotateFile.write(ANSI.FG_CYAN + subscript[last_end:m.end()] + ANSI.RESET)
                # prepend to result list, making ranges explicit:
                #    fortran has fastest-running index first
                #  while
                #    python/c has fastest-running index last
                result[0:0] = [ [int(m.group(3))] ]
                last_end = m.end()
                continue
            elif m.group(1) is not None:
                if self.annotateFile:
                    self.annotateFile.write(ANSI.FG_BRIGHT_CYAN + subscript[last_end:m.end()] + ANSI.RESET)
                # prepend to result list, making ranges explicit
                #    fortran has fastest-running index first
                #  while
                #    python/c has fastest-running index last
                result[0:0] = [list(range(int(m.group(1)),int(m.group(2))+1))]
                last_end = m.end()
                continue
            break
        if last_end < len(subscript):
            if subscript[last_end:].stript():
                LOGGER.error("ERROR: leftover chars in subscript: '%s'", subscript[last_end:])
                if self.annotateFile:
                    self.annotateFile.write(ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_RED + subscript[last_end:] + ANSI.RESET)
                self.bad_input = True
            else:
                if self.annotateFile:
                    self.annotateFile.write(ANSI.BEGIN_INVERT + ANSI.FG_BLUE + subscript[last_end:] + ANSI.RESET)
        return result

    def parse_line_state0(self, line, pos_in_line):
        # we have no open group
        m = cRE_start_group.match(line, pos_in_line)
        if m is not None:
            self.__nl_group = m.group(1).lower()
            if self.annotateFile:
                self.annotateFile.write(ANSI.FG_BRIGHT_GREEN + m.group() + ANSI.RESET)
            self.state = 1
            self.onOpen_namelist_group(self.__nl_group)
            return m.end()
        else:
            # but comments may appear here
            m = cRE_comment.match(line, pos_in_line)
            if m is not None:
                if self.annotateFile:
                    self.annotateFile.write(ANSI.FG_BLUE + m.group() + ANSI.RESET)
                self.onComment(m.group())
                return m.end()
        return None

    def parse_line_state3(self, line, pos_in_line):
        # we are inside single-quoted multiline string
        m = cRE_str_s_close.match(line, pos_in_line)
        if m is None:
            if self.annotateFile:
                self.annotateFile.write(ANSI.FG_YELLOW + line[pos_in_line:] + ANSI.RESET)
            self.__values[-1] += "\n" + line
            return len(line)
        else:
            if self.annotateFile:
                self.annotateFile.write(ANSI.FG_YELLOW + m.group() + ANSI.RESET)
            self.__values[-1] += "\n" + m.group(1)
            self.__values[-1] = unquote_string(self.__values[-1])
            self.__types[-1] = 'C'
            self.state = 2
            return m.end()
        return None

    def parse_line_state4(self, line, pos_in_line):
        # we are inside double-quoted multiline string
        m = cRE_str_d_close.match(line, pos_in_line)
        if m is None:
            if self.annotateFile:
                self.annotateFile.write(ANSI.FG_YELLOW + line[pos_in_line:] + ANSI.RESET)
            self.__values[-1] += "\n" + line
            return len(line)
        else:
            if self.annotateFile:
                self.annotateFile.write(ANSI.FG_YELLOW + m.group() + ANSI.RESET)
            self.__values[-1] += "\n" + m.group(1)
            self.__values[-1] = unquote_string(self.__values[-1])
            self.__types[-1] = 'C'
            self.state = 2
            return m.end()
        return None

    def parse_line(self, line):
        pos_in_line = 0
        while pos_in_line<len(line):
            new_pos_in_line = None
            if self.state == 0:
                new_pos_in_line = self.parse_line_state0(line, pos_in_line)
            elif self.state==3:
                new_pos_in_line = self.parse_line_state3(line, pos_in_line)
            elif self.state==4:
                new_pos_in_line = self.parse_line_state4(line, pos_in_line)
            elif self.state == 1 or self.state == 2:
                # we are inside opened group
                #   check for group-closing /
                m = cRE_end_group.match(line, pos_in_line)
                if m is not None:
                    if self.__target is not None:
                        self.onClose_value_assignment(
                            self.__nl_group,
                            self.__target, self.__subscript,
                            self.__values, self.__types)
                    self.__target = None
                    self.__subscript = None
                    self.__values = None
                    self.__types = None
                    self.__nvalues_after_comma = 0
                    self.onClose_namelist_group(self.__nl_group)
                    self.__nl_group = None
                    if self.annotateFile:
                        self.annotateFile.write(ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_GREEN + m.group() + ANSI.RESET)
                    self.state = 0
                    new_pos_in_line = m.end()
                #   check for new assignment
                m = cRE_start_assignment.match(line, pos_in_line)
                if m is not None:
                    if self.__target is not None:
                        self.onClose_value_assignment(
                            self.__nl_group,
                            self.__target, self.__subscript,
                            self.__values, self.__types)
                    self.state = 2
                    if m.group('subscript') is None:
                        self.__subscript = None
                        if self.annotateFile:
                            self.annotateFile.write(ANSI.FG_GREEN + m.group() + ANSI.RESET)
                    else:
                        if self.annotateFile:
                            self.annotateFile.write(ANSI.FG_GREEN + line[pos_in_line:m.start('subscript')] + ANSI.RESET)
                        self.__subscript = self.parse_subscript(m.group('subscript'))
                        if self.annotateFile:
                            self.annotateFile.write(ANSI.FG_GREEN + line[m.end('subscript'):m.end()] + ANSI.RESET)
                    self.__target = m.group('target').lower()
                    self.__values = []
                    self.__types = []
                    self.__nvalues_after_comma = 0
                    self.onOpen_value_assignment(
                        self.__nl_group,
                        self.__target, self.__subscript)
                    new_pos_in_line=m.end()
                if self.state == 2:
                    # we are inside the values-part of an assignment
                    m = cRE_assigned_value.match(line, pos_in_line)
                    if m is not None:
                        new_pos_in_line=m.end()
                        if m.group('comment') is not None:
                            if self.annotateFile:
                                self.annotateFile.write(ANSI.FG_BLUE + m.group() + ANSI.RESET)
                            self.onComment(m.group())
                        else:
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
                            self.__nvalues_after_comma +=1
                            if self.annotateFile:
                                self.annotateFile.write(ANSI.FG_YELLOW + m.group() + ANSI.RESET)
                    else:
                        # special meaning of comma: may indicate Null values in array assignments
                        m = cRE_comma.match(line, pos_in_line)
                        if m is not None:
                            if self.__nvalues_after_comma is 0:
                                self.__values.append(None)
                                self.__types.append(None)
                            self.__nvalues_after_comma = 0
                            if self.annotateFile:
                                self.annotateFile.write(ANSI.FG_MAGENTA + m.group() + ANSI.RESET)
                            new_pos_in_line = m.end()
            if new_pos_in_line is None:
                break
            else:
                pos_in_line = new_pos_in_line
        if pos_in_line < len(line):
            if self.state > 0 and self.state < 5 and line[pos_in_line:].strip():
                # states we as the base class are handling, but with leftover chars on a line
                LOGGER.error("ERROR: leftover chars in line while inside namelist group: '%s'", line[pos_in_line:])
                self.bad_input = True
                if self.annotateFile:
                    self.annotateFile.write(ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_RED + line[pos_in_line:] + ANSI.RESET)
            else:
                if self.annotateFile:
                    self.annotateFile.write(ANSI.BEGIN_INVERT + ANSI.FG_BLUE + line[pos_in_line:] + ANSI.RESET)
        if self.annotateFile:
            self.annotateFile.write('\n')

    # Hooks to be overloaded in derived classes in order to do stuff
    def onComment(self, comment):
        pass

    def onOpen_namelist_group(self, groupname):
        pass

    def onClose_namelist_group(self, groupname):
        pass

    def onOpen_value_assignment(self, groupname, target, subscript):
        pass

    def onClose_value_assignment(self, groupname, target, subscript, values, dtypes):
        if subscript is None:
            LOGGER.error("SET %s/%s = %s (types: %s)", groupname, target, str(values), str(dtypes))
        else:
            LOGGER.error("SET %s/%s(%s) = %s (types: %s)", groupname, target, subscript, str(values), str(dtypes))

    def onBad_input(self):
        pass

if __name__ == "__main__":
    parser = FortranNamelistParser(sys.argv[1])
    parser.parse()
