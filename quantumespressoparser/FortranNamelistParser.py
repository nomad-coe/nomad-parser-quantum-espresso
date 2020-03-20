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
    if dtype == 'f':
        return (float(pyfloat_str), dtype)
    else:
        return (int(pyfloat_str), dtype)

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
cRE_comment = re.compile(r"\s*!(?P<comment>.*)")
RE_identifier = r"[a-zA-Z]\w*" # fortran identifier
cRE_start_group = re.compile(r'\s*&(' + RE_identifier + r')') # beginning of namelist group
cRE_end_group = re.compile(r'\s*/')
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
        r'!(?P<comment>.*)', # comment
    ]) + ')', re.I)
cRE_str_s_close = re.compile(r"([^']*(?:[^']|'')*'(?!'))") # single-quoted string, closing
cRE_str_d_close = re.compile(r'([^"]*(?:[^"]|"")*"(?!"))') # double-quoted string, closing
cRE_comma = re.compile(r'\s*,')
cRE_trailing_whitespace = re.compile(r'\s+$')

cRE_identifier = re.compile(r'\s*(?P<target>' + RE_identifier + r')')
cRE_assignment_subscript_open = re.compile(r'\s*\((?P<subscript>[^\)!]*)')
cRE_assignment_subscript_continue = re.compile(r'(?P<subscript>[^\)!]+)')
cRE_assignment_subscript_close = re.compile(r'(?P<subscript>[^\)!]*)\)')
cRE_assignment_equals = re.compile(r'\s*=')

cRE_subscript = re.compile(r'\s*,?\s*(?:(\d*)\s*:\s*(\d*)|(\d+))')

cRE_end_newline = re.compile(r'(.*?)(\n*)$')

class FortranNamelistParser(object):
    """Parser for Fortran 90 Namelists
    """
    def __init__(self, file_path, annotateFile = None):
        self.input_tree = {}
        self.file_path = file_path
        self.state = self.state_root
        self.__annotateFile = annotateFile
        self.__nl_group = None
        self.__target = None
        self.__subscript = None
        self.__values = None
        self.__types = None
        self.__nvalues_after_comma = 0
        self.__cre_closing = None
        self.bad_input = False
        self.cache = {}

    def parse(self):
        """open file and parse line-by-line"""
        with open(self.file_path, "r") as fIn:
            # process line-by-line
            for line in fIn:
                self.parse_line(line)
        # check if there was input flagged as 'bad'/'syntactically incorrect'
        if self.bad_input:
            # call bad-input hook
            self.onBad_input()
        # call end-of-file hook
        self.onEnd_of_file()

    def parse_line(self, line):
        """parse one line, delegating to the parser state handlers"""
        pos_in_line = 0
        while pos_in_line<len(line):
            new_pos_in_line = self.state(line, pos_in_line)
            # check if anything was parsed, otherwise cancel that line
            if new_pos_in_line is None:
                break
            else:
                pos_in_line = new_pos_in_line
        if pos_in_line < len(line):
            self.bad_input = True
            self.annotate(line[pos_in_line:], ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_RED)

    def annotate(self, what, highlight):
        """write string to annotateFile with ANSI highlight/reset sequences"""
        if self.__annotateFile:
            m = cRE_end_newline.match(what)
            self.__annotateFile.write(highlight + m.group(1) + ANSI.RESET + m.group(2))

    def parse_subscript_string(self, subscript):
        """parse fully captured subscript string into python array"""
        if subscript is None:
            return None
        result = []
        last_end = 0
        while last_end<len(subscript):
            m = cRE_subscript.match(subscript, last_end)
            if m is None:
                break
            elif m.group(3) is not None:
                # prepend to result list, making ranges explicit:
                #    fortran has fastest-running index first
                #  while
                #    python/c has fastest-running index last
                result[0:0] = [ [int(m.group(3))] ]
                last_end = m.end()
                continue
            elif m.group(1) is not None:
                # prepend to result list, making ranges explicit
                #    fortran has fastest-running index first
                #  while
                #    python/c has fastest-running index last
                result[0:0] = [list(range(int(m.group(1)),int(m.group(2))+1))]
                last_end = m.end()
                continue
            break
        if last_end < len(subscript):
            if subscript[last_end:].strip():
                LOGGER.error("ERROR: leftover chars in subscript: '%s'", subscript[last_end:])
                self.bad_input = True
        return result

    def state_root(self, line, pos_in_line):
        """state: no open namelist groups, i.e. at the root of the namelist"""
        m = cRE_start_group.match(line, pos_in_line)
        if m is not None:
            self.__nl_group = m.group(1).lower()
            self.annotate(m.group(), ANSI.FG_BRIGHT_GREEN)
            self.state = self.state_inside_group
            self.onOpen_namelist_group(self.__nl_group)
            return m.end()
        else:
            # but comments may appear here
            m = cRE_comment.match(line, pos_in_line)
            if m is not None:
                self.annotate(m.group(), ANSI.FG_BLUE)
                self.onComment(m.group('comment'))
                return m.end()
            # as well as whitespace-only lines
            m = cRE_trailing_whitespace.match(line, pos_in_line)
            if m is not None:
                self.annotate(m.group(), ANSI.BG_WHITE)
                return m.end()
        # nothing matched, call hook
        return self.onRoot_data(line, pos_in_line)

    def state_inside_group(self, line, pos_in_line):
        """state: inside opened group, but no open assignment"""
        # check for group-closing /
        m = cRE_end_group.match(line, pos_in_line)
        if m is not None:
            # we just closed a NL group
            self.annotate(m.group(), ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_GREEN)
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
            self.state = self.state_root
            return m.end()
        # check for new identifier (part of left-hand side of assignment)
        m = cRE_identifier.match(line, pos_in_line)
        if m is not None:
            self.annotate(m.group(), ANSI.FG_GREEN)
            if self.__target is not None:
                self.onClose_value_assignment(
                    self.__nl_group,
                    self.__target, self.__subscript,
                    self.__values, self.__types)
            self.__target = m.group('target').lower()
            self.__subscript = None
            self.__values = []
            self.__types = []
            self.__nvalues_after_comma = 0
            return m.end()
        # check for new subscript (part of left-hand side of assignment)
        m = cRE_assignment_subscript_open.match(line, pos_in_line)
        if m is not None:
            self.annotate(line[pos_in_line:m.start('subscript')], ANSI.FG_GREEN)
            self.annotate(m.group('subscript'), ANSI.FG_CYAN)
            self.__subscript = m.group('subscript')
            self.state = self.state_assignment_subscript
            return m.end()
        # check for '=' sign in assignment, separating target from values
        m = cRE_assignment_equals.match(line, pos_in_line)
        if m is not None:
            self.annotate(line[pos_in_line:m.end()], ANSI.FG_GREEN)
            self.onOpen_value_assignment(
                self.__nl_group,
                self.__target, self.__subscript)
            self.state = self.state_assignment_values
            return m.end()
        # check for comments ('!' character up until end of line)
        m = cRE_comment.match(line, pos_in_line)
        if m is not None:
            self.annotate(m.group(), ANSI.FG_BLUE)
            self.onComment(m.group('comment'))
            return m.end()
        # check for trailing whitespace
        m = cRE_trailing_whitespace.match(line, pos_in_line)
        if m is not None:
            self.annotate(m.group(), ANSI.BG_WHITE)
            return m.end()
        return None

    def state_assignment_values(self, line, pos_in_line):
        """state: parse values, i.e. right-hand side of assignment"""
        # match value literals, groups decide on data type
        m = cRE_assigned_value.match(line, pos_in_line)
        if m is not None:
            if m.group('comment') is not None:
                # found a comment
                self.annotate(m.group(), ANSI.FG_BLUE)
                self.onComment(m.group('comment'))
            else:
                self.annotate(m.group(), ANSI.FG_YELLOW)
                if m.group('num') is not None:
                    # literal is single integer or float
                    (value, dtype) = match_to_float(m, group_offset=1)
                    self.__values.append(value)
                    self.__types.append(dtype)
                elif m.group('cnum_r') is not None:
                    # literal is complex: (float, float)
                    (cnum_r, dtype) = match_to_float(m, group_offset=10)
                    (cnum_i, dtype) = match_to_float(m, group_offset=19)
                    self.__values.append(complex(cnum_r, cnum_i))
                    self.__types.append('complex')
                elif m.group('bool_t') is not None:
                    # literal is a true-value bool
                    self.__values.append(True)
                    self.__types.append('b')
                elif m.group('bool_f') is not None:
                    # literal is a false-value bool
                    self.__values.append(False)
                    self.__types.append('b')
                elif m.group('str_s') is not None:
                    # literal is a closed, single-quoted string
                    self.__values.append(unquote_string(m.group('str_s')))
                    self.__types.append('C')
                elif m.group('str_d') is not None:
                    # literal is a closed, double-quoted string
                    self.__values.append(unquote_string(m.group('str_d')))
                    self.__types.append('C')
                elif m.group('str_s_nc') is not None:
                    # literal is a non-closed, single-quoted string
                    self.state = self.state_assignment_values_multiline_string
                    self.__values.append(m.group('str_s_nc'))
                    self.__types.append('C')
                    self.__cre_closing = cRE_str_s_close
                elif m.group('str_d_nc') is not None:
                    # literal is a non-closed, double-quoted string
                    self.state = self.state_assignment_values_multiline_string
                    self.__values.append(m.group('str_d_nc'))
                    self.__types.append('C')
                    self.__cre_closing = cRE_str_d_close
                # keep track if there were values following the previous comma
                self.__nvalues_after_comma += 1
            return m.end()
        # special meaning of comma: may indicate Null values in array assignments
        m = cRE_comma.match(line, pos_in_line)
        if m is not None:
            self.annotate(m.group(), ANSI.FG_MAGENTA)
            if self.__nvalues_after_comma is 0:
                # there were no value literals preceeding the comma, which
                # means a null value
                self.__values.append(None)
                self.__types.append(None)
            self.__nvalues_after_comma = 0
            return m.end()
        # check for trailing whitespace
        m = cRE_trailing_whitespace.match(line, pos_in_line)
        if m is not None:
            self.annotate(m.group(), ANSI.BG_WHITE)
            return m.end()
        # if none of the above matched, switch back to checking for new assignment
        self.state = self.state_inside_group
        return pos_in_line

    def state_assignment_values_multiline_string(self, line, pos_in_line):
        """state: parse multiline string in right-hand side of assignment"""
        # check for closing quotes
        m = self.__cre_closing.match(line, pos_in_line)
        if m is None:
            # no closing quotes, append data to value
            self.annotate(line[pos_in_line:], ANSI.FG_YELLOW)
            self.__values[-1] += line
            return len(line)
        else:
            # closing quotes, postprocess string
            self.annotate(m.group(), ANSI.FG_YELLOW)
            self.__values[-1] += m.group(1)
            # remove enclosing quotes and resolve escaped quotes in string
            self.__values[-1] = unquote_string(self.__values[-1])
            self.__cre_closing = None
            self.state = self.state_assignment_values
            return m.end()
        return None

    def state_assignment_subscript(self, line, pos_in_line):
        """state: capture subscipt, possibly spanning multiple lines"""
        # check for closing bracket
        m = cRE_assignment_subscript_close.match(line, pos_in_line)
        if m is not None:
            # subscript closed, convert string form to python array
            self.annotate(m.group('subscript'), ANSI.FG_CYAN)
            self.annotate(line[m.end('subscript'):m.end()], ANSI.FG_GREEN)
            self.__subscript = self.parse_subscript_string(self.__subscript + m.group('subscript'))
            self.state = self.state_inside_group
            return m.end()
        # check for new indices in subscript
        m = cRE_assignment_subscript_continue.match(line, pos_in_line)
        if m is not None:
            self.annotate(m.group('subscript'), ANSI.FG_CYAN)
            self.__subscript += m.group('subscript')
            return m.end()
        # comments may appear within subscripts spanning multiple lines
        m = cRE_comment.match(line, pos_in_line)
        if m is not None:
            self.annotate(m.group(), ANSI.FG_BLUE)
            self.onComment(m.group('comment'))
            self.__subscript += line[pos_in_line:m.start()]
            return m.end()
        self.annotate(line[pos_in_line:], ANSI.BEGIN_INVERT + ANSI.FG_BRIGHT_RED)
        LOGGER.error("ERROR: leftover chars in line while inside subscript: '%s'", line[pos_in_line:])
        self.bad_input = True
        return None

    # Hooks to be overloaded in derived classes in order to do stuff beyond caching
    def onComment(self, comment):
        """hook: called whan a comment was found"""
        pass

    def onOpen_namelist_group(self, groupname):
        """hook: called when a namelist group opens"""
        if groupname in self.cache:
            LOGGER.error("ERROR: multiple definitions of group &%s", groupname)
            self.bad_input = True
        else:
            self.cache[groupname]={}

    def onClose_namelist_group(self, groupname):
        """hook: called when a namelist group closes"""
        LOGGER.error("group: %s", groupname)
        for identifier in sorted(self.cache[groupname]):
            LOGGER.error("  %s: %s", identifier, str(self.cache[groupname][identifier]))

    def onOpen_value_assignment(self, groupname, target, subscript):
        """hook: called when a value assignment within a namelist group starts"""
        pass

    def onClose_value_assignment(self, groupname, target, subscript, values, dtypes):
        """hook: called when a value assignment within a namelist group closes
        Arguments are: NL group name, identifier/subscript, values and assumed
        data types
        """
        if subscript is None:
            LOGGER.debug("SET %s/%s = %s (types: %s)", groupname, target, str(values), str(dtypes))
        else:
            LOGGER.debug("SET %s/%s(%s) = %s (types: %s)", groupname, target, subscript, str(values), str(dtypes))
        if target not in self.cache[groupname]:
            self.cache[groupname][target] = []
        self.cache[groupname][target].append([subscript, values, dtypes])

    def onRoot_data(self, line, pos_in_line):
        """hook: called if data appears outside namelists groups, directly
        at root level within the file;
        data means: line is not empty or a comment
        useful for code-specific extensions beyond the F90 namelist standard
        """
        return None

    def onBad_input(self):
        """hook: called at the end of parsing if there was any bad input"""
        pass

    def onEnd_of_file(self):
        """hook: called at the end of parsing"""
        pass

if __name__ == "__main__":
    parser = FortranNamelistParser(sys.argv[1], annotateFile=sys.stdout)
    parser.parse()
