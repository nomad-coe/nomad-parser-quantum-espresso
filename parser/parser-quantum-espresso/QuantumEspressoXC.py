import logging
import re
import copy


LOGGER = logging.getLogger(__name__)


def parse_qe_xc_num(xc_functional_num):
    """parse Quantum Espresso XC number/index notation
    :returns: list with 6 integer elements (components of QE XC functionals)
              [iexch, icorr, igcx, igcc, imeta, inlc]
                iexch, icorr - density exchange/correlation
                igcx, igcc   - gradient correction exchange/correlation
                imeta        - meta-GGA
                inlc         - non-local term of Van der Waals functionals
    """
    xf_num_split_i = []
    xf_num_split = re.split(r'\s+', xc_functional_num.strip())
    if len(xf_num_split) > 6:
        # tested versions of espresso <=5.4.0 have 6 elements
        raise RuntimeError("unsupported number of XC components: %d" % (
            len(xf_num_split)))
    elif len(xf_num_split) == 6:
        # trivial case: we got 6 components from simply splitting
        xf_num_split_i = [int(x) for x in xf_num_split]
    elif len(xf_num_split) == 1 and len(xc_functional_num) == 4:
        # old 4-component form without space separator
        xf_num_split_i = [int(x) for x in re.findall(
            '(\d)', xc_functional_num)]
    elif len(xc_functional_num) == 10:
        # intermediate versions used 2-digit, 5 component form,
        # occasionally missing spaces
        # examples:
        # ( 1 4 4 0 1)
        # ( 0 413 0 2)  <- missing space
        xf_num_split_i = [int(x) for x in re.findall(
            '([ \d]\d)', xc_functional_num)]
    else:
        raise RuntimeError("unparsable input: '%s'", xc_functional_num)
    if len(xf_num_split_i) < 1:
        raise RuntimeError("this should not happen")
    # zero-pad up to 6 elements
    xf_num_split_i += [0] * (6 - len(xf_num_split_i))
    return xf_num_split_i


def translate_qe_xc_num(xc_functional_num, exact_exchange_fraction=None):
    # all remainders of excact_exchange_fraction must be ordinary DFT
    if exact_exchange_fraction is None:
        exact_exchange_fraction = 0.0
        dft_exchange_fraction = 1.0
    else:
        dft_exchange_fraction = 1.0 - exact_exchange_fraction
    xf_num = parse_qe_xc_num(xc_functional_num)
    LOGGER.debug('num <- input: %s <- %s, exx_fraction: %s',  str(xf_num),
                 xc_functional_num, str(exact_exchange_fraction))
    # use dictionary to ensure uniqueness:
    #   exchange/correlation functionals may be combined into _XC_, and we
    #   only want to emit such combinations once
    xc_data = {}
    xc_data_subtract = {}
    xc_section_method = {}
    for component_i in range(6):
        this_xf_num = xf_num[component_i]
        if this_xf_num == 0:
            # 0 means unset component
            continue
        component_max = len(XC_COMPONENT[component_i])-1
        this_component = None
        if this_xf_num > component_max:
            LOGGER.error(
                "%s[%d] beyond limit of %d",
                XC_COMPONENT_NAME[component_i], this_xf_num, component_max)
        else:
            this_component = XC_COMPONENT[component_i][this_xf_num]
        if this_component is None:
            raise RuntimeError("Undefined XC component %s[%d]" % (
                XC_COMPONENT_NAME[component_i], this_xf_num))
        xc_section_method.update(this_component['xc_section_method'])
        for this_term in this_component['xc_terms']:
            apply_term_add(
                xc_data, this_term,
                exact_exchange_fraction, dft_exchange_fraction)
        if 'xc_terms_subtract' in this_component:
            for this_term in this_component['xc_terms_subtract']:
                apply_term_add(
                    xc_data_subtract, this_term,
                    exact_exchange_fraction, dft_exchange_fraction)
    apply_terms_subtract(xc_data, xc_data_subtract)
    apply_terms_filter(xc_data)
    result = []
    for (k, v) in sorted(xc_data.items()):
        result.append(v)
    return (xc_section_method, result)


def apply_term_add(xc_data, this_term,
                   exact_exchange_fraction, dft_exchange_fraction):
    term = copy.deepcopy(this_term)
    if term['XC_functional_name'] == 'HYB_GGA_XC_HSE06':
        if exact_exchange_fraction is not None:
            # we are at HSE06 component, with explicit exact_exchange_fraction
            term['XC_functional_parameters'] = {
                'exx_mixing': exact_exchange_fraction,
            }
    if term.get('t_qe_XC_functional_weight_scale_exx', None):
        term['XC_functional_weight'] = (
            term['t_qe_XC_functional_weight_scale_exx'] *
            exact_exchange_fraction
        )
    if term.get('t_qe_XC_functional_weight_scale_dft', None):
        term['XC_functional_weight'] = (
            term['t_qe_XC_functional_weight_scale_dft'] *
            dft_exchange_fraction
        )
    if 'XC_functional_weight' not in term:
        term['XC_functional_weight'] = 1.0
    if term['XC_functional_name'] not in xc_data:
        xc_data[term['XC_functional_name']] = term
    else:
        LOGGER.info("pre-existing XC term: %s",
                    term['XC_functional_name'])
    return xc_data


def apply_terms_subtract(xc_data, xc_data_remove):
    for (k, v) in xc_data_remove.items():
        if k in xc_data:
            xc_data[k]['XC_functional_weight'] -= v['XC_functional_weight']
        else:
            xc_data[k] = v
            xc_data[k]['XC_functional_weight'] *= -1.0


def apply_terms_filter(xc_data):
    for (k, v) in list(xc_data.items()):
        if abs(v['XC_functional_weight']) < 0.01:
            del xc_data[k]
        else:
            if abs(v['XC_functional_weight'] - 1.0) < 0.01:
                del v['XC_functional_weight']
            v.pop('t_qe_XC_functional_weight_scale_exx', None)
            v.pop('t_qe_XC_functional_weight_scale_dft', None)


# origin: espresso-5.4.0/Modules/funct.f90
EXCHANGE = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       'sla',
            'x_qe_xc_iexch_comment':    'Slater (alpha=2/3)',
            'x_qe_xc_iexch':      1,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_X',
            'XC_functional_parameters': {'alpha': 1.0},
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       'sl1',
            'x_qe_xc_iexch_comment':    'Slater (alpha=1.0)',
            'x_qe_xc_iexch':      2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'x_qe_LDA_X_RELATIVISTIC',
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       'rxc',
            'x_qe_xc_iexch_comment':    'Relativistic Slater',
            'x_qe_xc_iexch':      3,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'OEP_EXX',
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       'oep',
            'x_qe_xc_iexch_comment':    'Optimized Effective Potential',
            'x_qe_xc_iexch':      4,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'HF_X',
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       'hf',
            'x_qe_xc_iexch_comment':    'Hartree-Fock',
            'x_qe_xc_iexch':      5,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'HF_X',
            't_qe_XC_functional_weight_scale_exx': 1.0,
            'XC_functional_weight': 0.25,
        }, {
            'XC_functional_name': 'LDA_X',
            't_qe_XC_functional_weight_scale_dft': 1.0,
            'XC_functional_weight': 0.75,
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       "pb0x",
            'x_qe_xc_iexch_comment':    'PBE0 (Slater*0.75+HF*0.25)',
            'x_qe_xc_iexch':      6,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'HF_X',
            't_qe_XC_functional_weight_scale_exx': 1.0,
            'XC_functional_weight': 0.20,
        }, {
            'XC_functional_name': 'LDA_X',
            't_qe_XC_functional_weight_scale_dft': 1.0,
            'XC_functional_weight': 0.8,
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       "b3lp",
            'x_qe_xc_iexch_comment':    "B3LYP(Slater*0.80+HF*0.20)",
            'x_qe_xc_iexch':      7,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_LDA_X_KZK",
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       "kzk",
            'x_qe_xc_iexch_comment':    "Finite-size corrections",
            'x_qe_xc_iexch':      8,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'HF_X',
            't_qe_XC_functional_weight_scale_exx': 1.0,
            'XC_functional_weight': 0.218,
        }, {
            'XC_functional_name': 'LDA_X',
            't_qe_XC_functional_weight_scale_dft': 1.0,
            'XC_functional_weight': 0.782,
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       "x3lp",
            'x_qe_xc_iexch_comment':    "X3LYP(Slater*0.782+HF*0.218)",
            'x_qe_xc_iexch':      9,
        },
    },
]


CORRELATION = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_PZ',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "pz",
            'x_qe_xc_icorr_comment':    "Perdew-Zunger",
            'x_qe_xc_icorr':      1,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_VWN',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "vwn",
            'x_qe_xc_icorr_comment':    "Vosko-Wilk-Nusair",
            'x_qe_xc_icorr':      2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_LYP',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "lyp",
            'x_qe_xc_icorr_comment':    "Lee-Yang-Parr",
            'x_qe_xc_icorr':      3,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "LDA_C_PW",
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "pw",
            'x_qe_xc_icorr_comment':    "Perdew-Wang",
            'x_qe_xc_icorr':      4,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_WIGNER',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "wig",
            'x_qe_xc_icorr_comment':    "Wigner",
            'x_qe_xc_icorr':      5,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_HL',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "hl",
            'x_qe_xc_icorr_comment':    "Hedin-Lunqvist",
            'x_qe_xc_icorr':      6,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_OB_PZ',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "obz",
            'x_qe_xc_icorr_comment':    "Ortiz-Ballone form for PZ",
            'x_qe_xc_icorr':      7,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_OB_PW',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "obw",
            'x_qe_xc_icorr_comment':    "Ortiz-Ballone form for PW",
            'x_qe_xc_icorr':      8,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_GL',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "gl",
            'x_qe_xc_icorr_comment':    "Gunnarson-Lunqvist",
            'x_qe_xc_icorr':      9,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_LDA_C_KZK",
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "kzk",
            'x_qe_xc_icorr_comment':    "Finite-size corrections",
            'x_qe_xc_icorr':      10,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_VWN_RPA',
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "vwn-rpa",
            'x_qe_xc_icorr_comment':    "Vosko-Wilk-Nusair, alt param",
            'x_qe_xc_icorr':      11,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_VWN',
            'XC_functional_weight': 0.19,
        }, {
            'XC_functional_name': 'LDA_C_LYP',
            'XC_functional_weight': 0.81,
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "b3lp",
            'x_qe_xc_icorr_comment':    "B3LYP (0.19*vwn+0.81*lyp)",
            'x_qe_xc_icorr':      12,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_LDA_C_VWN_RPA',
            'XC_functional_weight': 0.19,
        }, {
            'XC_functional_name': 'LDA_C_LYP',
            'XC_functional_weight': 0.81,
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':    "b3lpv1r",
            'x_qe_xc_icorr_comment': "B3LYP-VWN-1-RPA (0.19*vwn_rpa+0.81*lyp)",
            'x_qe_xc_icorr':         13,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': 'LDA_C_LDA_C_VWN_RPA',
            'XC_functional_weight': 0.129,
        }, {
            'XC_functional_name': 'LDA_C_LYP',
            'XC_functional_weight': 0.871,
        }],
        'xc_section_method': {
            'x_qe_xc_icorr_name':       "x3lp",
            'x_qe_xc_icorr_comment':    "X3LYP (0.129*vwn_rpa+0.871*lyp)",
            'x_qe_xc_icorr':      14,
        },
    },
]

EXCHANGE_GRADIENT_CORRECTION = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B88",
        }],
        'xc_terms_subtract': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "b88",
            'x_qe_xc_igcx_comment':    "Becke88 (beta=0.0042)",
            'x_qe_xc_igcx':      1,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PW91",
        }],
        'xc_terms_subtract': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "ggx",
            'x_qe_xc_igcx_comment':    "Perdew-Wang 91",
            'x_qe_xc_igcx':      2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PBE",
        }],
        'xc_terms_subtract': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "pbx",
            'x_qe_xc_igcx_comment':    "Perdew-Burke-Ernzenhof exch",
            'x_qe_xc_igcx':      3,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PBE_R",
        }],
        'xc_terms_subtract': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "rpb",
            'x_qe_xc_igcx_comment':    "revised PBE by Zhang-Yang",
            'x_qe_xc_igcx':      4,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_XC_HCTH_120",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "hcth",
            'x_qe_xc_igcx_comment':    "Cambridge exch, Handy et al, HCTH/120",
            'x_qe_xc_igcx':      5,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_OPTX",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "optx",
            'x_qe_xc_igcx_comment':    "Handy's exchange functional",
            'x_qe_xc_igcx':      6,
        },
    },
    {
        # igcx=7 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_TPSS",
        }],
        'xc_terms_subtract': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':    "tpss",
            'x_qe_xc_igcx_comment': "TPSS Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcx':         7,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PBE",
            'XC_functional_weight': 0.75,
            't_qe_XC_functional_weight_scale_dft': 1.0,
        }],
        'xc_terms_subtract': [{
            'XC_functional_name': 'LDA_X',
            'XC_functional_weight': 0.75,
            't_qe_XC_functional_weight_scale_dft': 1.0,
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "pb0x",
            'x_qe_xc_igcx_comment':    "PBE0 (PBE exchange*0.75)",
            'x_qe_xc_igcx':      8,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B88",
            'XC_functional_weight': 0.72,
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "b3lp",
            'x_qe_xc_igcx_comment':    "B3LYP (Becke88*0.72)",
            'x_qe_xc_igcx':      9,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PBE_SOL",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "psx",
            'x_qe_xc_igcx_comment':    "PBEsol exchange",
            'x_qe_xc_igcx':      10,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_WC",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "wcx",
            'x_qe_xc_igcx_comment':    "Wu-Cohen",
            'x_qe_xc_igcx':      11,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "HYB_GGA_XC_HSE06",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "hse",
            'x_qe_xc_igcx_comment':    "HSE screened exchange",
            'x_qe_xc_igcx':      12,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RPW86",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "rw86",
            'x_qe_xc_igcx_comment':    "revised PW86",
            'x_qe_xc_igcx':      13,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PBE",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "pbe",
            'x_qe_xc_igcx_comment':    "same as PBX, back-comp.",
            'x_qe_xc_igcx':      14,
        },
    },
    {
        # igcx=15 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_TB09",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':    "tb09",
            'x_qe_xc_igcx_comment': "TB09 Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcx':         15,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_C09X",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "c09x",
            'x_qe_xc_igcx_comment':    "Cooper 09",
            'x_qe_xc_igcx':      16,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_SOGGA",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "sox",
            'x_qe_xc_igcx_comment':    "sogga",
            'x_qe_xc_igcx':      17,
        },
    },
    {
        # igcx=18 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_M06_L",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':    "m6lx",
            'x_qe_xc_igcx_comment': "M06L Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcx':         18,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_Q2D",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "q2dx",
            'x_qe_xc_igcx_comment':    "Q2D exchange grad corr",
            'x_qe_xc_igcx':      19,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_GGA_X_GAUP",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "gaup",
            'x_qe_xc_igcx_comment':    "Gau-PBE hybrid exchange",
            'x_qe_xc_igcx':      20,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PW86",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "pw86",
            'x_qe_xc_igcx_comment':    "Perdew-Wang (1986) exchange",
            'x_qe_xc_igcx':      21,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B86_MGC",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "b86b",
            'x_qe_xc_igcx_comment':    "Becke (1986) exchange",
            'x_qe_xc_igcx':      22,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_OPTB88_VDW",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "obk8",
            'x_qe_xc_igcx_comment':    "optB88  exchange",
            'x_qe_xc_igcx':      23,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_GGA_X_OPTB86_VDW",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "ob86",
            'x_qe_xc_igcx_comment':    "optB86b exchange",
            'x_qe_xc_igcx':      24,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_EV93",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "evx",
            'x_qe_xc_igcx_comment':    "Engel-Vosko exchange",
            'x_qe_xc_igcx':      25,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B86_R",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "b86r",
            'x_qe_xc_igcx_comment':    "revised Becke (b86b)",
            'x_qe_xc_igcx':      26,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_LV_RPW86",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "cx13",
            'x_qe_xc_igcx_comment':    "consistent exchange",
            'x_qe_xc_igcx':      27,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B88",
            'XC_functional_weight': 0.542,
        }, {
            'XC_functional_name': "GGA_X_PW91",
            'XC_functional_weight': 0.167,
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':    "x3lp",
            'x_qe_xc_igcx_comment': "X3LYP (Becke88*0.542 "
                                    " + Perdew-Wang91*0.167)",
            'x_qe_xc_igcx':         28,
        },
    },
]

CORRELATION_GRADIENT_CORRECTION = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_P86",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "p86",
            'x_qe_xc_igcc_comment':    "Perdew86",
            'x_qe_xc_igcc':      1,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_PW91",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "ggc",
            'x_qe_xc_igcc_comment':    "Perdew-Wang 91 corr.",
            'x_qe_xc_igcc':      2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_LYP",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "blyp",
            'x_qe_xc_igcc_comment':    "Lee-Yang-Parr",
            'x_qe_xc_igcc':      3,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_PBE",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "pbc",
            'x_qe_xc_igcc_comment':    "Perdew-Burke-Ernzenhof corr",
            'x_qe_xc_igcc':      4,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_XC_HCTH_120",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "hcth",
            'x_qe_xc_igcc_comment':    "Cambridge exch, Handy et al, HCTH/120",
            'x_qe_xc_igcc':      5,
        },
    },
    {
        # igcc=6 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [{
            'XC_functional_name': "MGGA_C_TPSS",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':    "tpss",
            'x_qe_xc_igcc_comment': "TPSS Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcc':         6,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_LYP",
            'XC_functional_weight': 0.81,
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "b3lp",
            'x_qe_xc_igcc_comment':    "B3LYP (Lee-Yang-Parr*0.81)",
            'x_qe_xc_igcc':      7,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_PBE_SOL",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "psc",
            'x_qe_xc_igcc_comment':    "PBEsol corr",
            'x_qe_xc_igcc':      8,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_PBE",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "pbe",
            'x_qe_xc_igcc_comment':    "same as PBX, back-comp.",
            'x_qe_xc_igcc':      9,
        },
    },
    {
        # igcc=10 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        #        functionals.f90 tells that correlation is taken from tpss
        'xc_terms': [{
            'XC_functional_name': "MGGA_C_TPSS",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':    "tb09",
            'x_qe_xc_igcc_comment': "TB09 Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcc':         10,
        },
    },
    {
        # igcc=11 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [{
            'XC_functional_name': "MGGA_C_M06_L",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':    "m6lx",
            'x_qe_xc_igcc_comment': "M06L Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcc':         2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_Q2D",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "q2dc",
            'x_qe_xc_igcc_comment':    "Q2D correlation grad corr",
            'x_qe_xc_igcc':      12,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_LYP",
            'XC_functional_weight': 0.871,
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "x3lp",
            'x_qe_xc_igcc_comment':    "X3LYP (Lee-Yang-Parr*0.871)",
            'x_qe_xc_igcc':      13,
        },
    },
]


META_GGA = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_TPSS",
        }, {
            'XC_functional_name': "MGGA_C_TPSS",
        }],
        'xc_terms_subtract': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_imeta_name':       "tpss",
            'x_qe_xc_imeta_comment':    "TPSS Meta-GGA",
            'x_qe_xc_imeta':      1,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_M06_L",
        }, {
            'XC_functional_name': "MGGA_C_M06_L",
        }],
        'xc_section_method': {
            'x_qe_xc_imeta_name':       "m6lx",
            'x_qe_xc_imeta_comment':    "M06L Meta-GGA",
            'x_qe_xc_imeta':      2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_TB09",
        }, {
            # confirmed by looking into functionals.f90
            'XC_functional_name': "MGGA_C_TPSS",
        }],
        'xc_section_method': {
            'x_qe_xc_imeta_name':       "tb09",
            'x_qe_xc_imeta_comment':    "TB09 Meta-GGA",
            'x_qe_xc_imeta':      3,
        },
    },
]

VAN_DER_WAALS = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_DF1",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vdw1",
            'x_qe_xc_inlc_comment':    "vdW-DF1",
            'x_qe_xc_inlc':      1,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_DF2",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vdw2",
            'x_qe_xc_inlc_comment':    "vdW-DF2",
            'x_qe_xc_inlc':      2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_RVV10",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vv10",
            'x_qe_xc_inlc_comment':    "rVV10",
            'x_qe_xc_inlc':      3,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_DFX",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':    "vdwx",
            'x_qe_xc_inlc_comment': "vdW-DF-x (reserved Thonhauser,"
                                    " not implemented)",
            'x_qe_xc_inlc':         4,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_DFY",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':    "vdwy",
            'x_qe_xc_inlc_comment': "vdW-DF-y (reserved Thonhauser,"
                                    " not implemented)",
            'x_qe_xc_inlc':         5,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_DFZ",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':    "vdwz",
            'x_qe_xc_inlc_comment': "vdW-DF-z (reserved Thonhauser,"
                                    " not implemented)",
            'x_qe_xc_inlc':         6,
        },
    },
]


XC_COMPONENT = [
    EXCHANGE,
    CORRELATION,
    EXCHANGE_GRADIENT_CORRECTION,
    CORRELATION_GRADIENT_CORRECTION,
    VAN_DER_WAALS,
    META_GGA,
]


XC_COMPONENT_NAME = [
    'EXCHANGE',
    'CORRELATION',
    'EXCHANGE_GRADIENT_CORRECTION',
    'CORRELATION_GRADIENT_CORRECTION',
    'VAN_DER_WAALS',
    'META_GGA',
]
