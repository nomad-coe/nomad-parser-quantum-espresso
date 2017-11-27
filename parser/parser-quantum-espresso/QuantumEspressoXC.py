import logging
import re


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
    xf_num_split = re.split(r'\s+',xc_functional_num.strip())
    if len(xf_num_split) > 6:
        # tested versions of espresso <=5.4.0 have 6 elements
        raise RuntimeError("unsupported number of XC components: %d" % (len(xf_num_split)))
    elif len(xf_num_split) == 6:
        # trivial case: we got 6 components from simply splitting
        xf_num_split_i = [ int(x) for x in xf_num_split ]
    elif len(xf_num_split) == 1 and len(xc_functional_num)==4:
        # old 4-component form without space separator
        xf_num_split_i = [ int(x) for x in re.findall('(\d)', xc_functional_num) ]
    elif len(xc_functional_num)==10:
        # intermediate versions used 2-digit, 5 component form, occasionally missing spaces
        # examples:
        # ( 1 4 4 0 1)
        # ( 0 413 0 2)  <- missing space
        xf_num_split_i = [ int(x) for x in re.findall('([ \d]\d)',xc_functional_num) ]
    else:
        raise RuntimeError("unparsable input: '%s'", xc_functional_num)
    if len(xf_num_split_i)<1:
        raise RuntimeError("this should not happen")
    # zero-pad up to 6 elements
    xf_num_split_i += [ 0 ] * (6 - len(xf_num_split_i))
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
            this_component=XC_COMPONENT[component_i][this_xf_num]
        if this_component is None:
            raise RuntimeError("Undefined XC component %s[%d]" % (
                XC_COMPONENT_NAME[component_i], this_xf_num))
        xc_section_method.update(this_component['xc_section_method'])
        for this_term in this_component['xc_terms']:
            add_term(xc_data, this_term,
                     exact_exchange_fraction, dft_exchange_fraction)
    result = []
    for k,v in sorted(xc_data.items()):
        result.append(v)
    return (xc_section_method, result)


def add_term(xc_data, this_term,
             exact_exchange_fraction, dft_exchange_fraction):
    term = this_term.copy()
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
        del term['t_qe_XC_functional_weight_scale_exx']
    if term.get('t_qe_XC_functional_weight_scale_dft', None):
        term['XC_functional_weight'] = (
            term['t_qe_XC_functional_weight_scale_dft'] *
            dft_exchange_fraction
        )
        del term['t_qe_XC_functional_weight_scale_dft']
    if term['XC_functional_name'] not in xc_data:
        xc_data[term['XC_functional_name']] = term
    else:
        LOGGER.info("pre-existing XC term: %s",
                    term['XC_functional_name'])
    return xc_data


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
            'XC_functional_parameters': { 'alpha': 1.0 },
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
            'x_qe_xc_icorr_name':       "b3lpv1r",
            'x_qe_xc_icorr_comment':    "B3LYP-VWN-1-RPA (0.19*vwn_rpa+0.81*lyp)",
            'x_qe_xc_icorr':      13,
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
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "pbx",
            'x_qe_xc_igcx_comment':    "Perdew-Burke-Ernzenhof exch",
            'x_qe_xc_igcx':      3,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PBE_R'",
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
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "tpss",
            'x_qe_xc_igcx_comment':    "TPSS Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcx':      7,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_PBE",
            'XC_functional_weight': 0.75,
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
            'x_qe_xc_igcx_name':       "tb09",
            'x_qe_xc_igcx_comment':    "TB09 Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcx':      15,
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
            'x_qe_xc_igcx_name':       "m6lx",
            'x_qe_xc_igcx_comment':    "M06L Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcx':      18,
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
            'x_qe_xc_igcx_name':       "x3lp",
            'x_qe_xc_igcx_comment':    "X3LYP (Becke88*0.542 + Perdew-Wang91*0.167)",
            'x_qe_xc_igcx':      28,
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
            'x_qe_xc_igcc_name':       "tpss",
            'x_qe_xc_igcc_comment':    "TPSS Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcc':      6,
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
            'x_qe_xc_igcc_name':       "tb09",
            'x_qe_xc_igcc_comment':    "TB09 Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcc':      10,
        },
    },
    {
        # igcc=11 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [{
            'XC_functional_name': "MGGA_C_M06_L",
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "m6lx",
            'x_qe_xc_igcc_comment':    "M06L Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcc':      2,
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


META_GGA= [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_TPSS",
        }, {
            'XC_functional_name': "MGGA_C_TPSS",
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
            'x_qe_xc_inlc_name':       "vdwx",
            'x_qe_xc_inlc_comment':    "vdW-DF-x (reserved Thonhauser, not implemented)",
            'x_qe_xc_inlc':      4,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_DFY",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vdwy",
            'x_qe_xc_inlc_comment':    "vdW-DF-y (reserved Thonhauser, not implemented)",
            'x_qe_xc_inlc':      5,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "x_qe_VDW_DFZ",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vdwz",
            'x_qe_xc_inlc_comment':    "vdW-DF-z (reserved Thonhauser, not implemented)",
            'x_qe_xc_inlc':      6,
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


# origin: espresso-5.4.0/Modules/funct.f90
#  "dft" is the exchange-correlation functional label, described either
#  by short names listed below, or by a series of keywords (everything
#  is case-insensitive). "dft_shortname" contains one of the short names
#  listed below (deduced from from "dft" as read from input or PP files)
#
#  short name       complete name       Short description
#     "pz"    = "sla+pz"            = Perdew-Zunger LDA
#     "bp"    = "b88+p86"           = Becke-Perdew grad.corr.
#     "pw91"  = "sla+pw+ggx+ggc"    = PW91 (aka GGA)
#     "blyp"  = "sla+b88+lyp+blyp"  = BLYP
#     "pbe"   = "sla+pw+pbx+pbc"    = PBE
#     "revpbe"= "sla+pw+rpb+pbc"    = revPBE (Zhang-Yang)
#     "pw86pbe" = "sla+pw+pw86+pbc" = PW86 exchange + PBE correlation
#     "b86bpbe" = "sla+pw+b86b+pbc" = B86b exchange + PBE correlation
#     "pbesol"= "sla+pw+psx+psc"    = PBEsol
#     "q2d"   = "sla+pw+q2dx+q2dc"  = PBEQ2D
#     "hcth"  = "nox+noc+hcth+hcth" = HCTH/120
#     "olyp"  = "nox+lyp+optx+blyp" = OLYP
#     "wc"    = "sla+pw+wcx+pbc"    = Wu-Cohen
#     "sogga" = "sla+pw+sox+pbec"   = SOGGA
#     "optbk88"="sla+pw+obk8+p86"   = optB88
#     "optb86b"="sla+pw+ob86+p86"   = optB86
#     "ev93"  = "sla+pw+evx+nogc"   = Engel-Vosko
#     "tpss"  = "sla+pw+tpss+tpss"  = TPSS Meta-GGA
#     "m06l"  = "nox+noc+m6lx+m6lc" = M06L Meta-GGA
#     "tb09"  = "sla+pw+tb09+tb09"  = TB09 Meta-GGA
#     "pbe0"  = "pb0x+pw+pb0x+pbc"  = PBE0
#     "hse"   = "sla+pw+hse+pbc"    = Heyd-Scuseria-Ernzerhof (HSE 06, see note below)
#     "b3lyp" = "b3lp+b3lp+b3lp+b3lp"= B3LYP
#     "b3lypv1r"    = "b3lp+b3lpv1r+b3lp+b3lp"= B3LYP-VWN1-RPA
#     "x3lyp" = "x3lp+x3lp+x3lp+x3lp"= X3LYP
#     "vwn-rpa"     = "sla+vwn-rpa" = VWN LDA using vwn1-rpa parametriz
#     "gaupbe"= "sla+pw+gaup+pbc"   = Gau-PBE (also "gaup")
#     "vdw-df"       ="sla+pw+rpb +vdw1"   = vdW-DF1
#     "vdw-df2"      ="sla+pw+rw86+vdw2"   = vdW-DF2
#     "vdw-df-x"     ="sla+pw+????+vdwx"   = vdW-DF-x, reserved Thonhauser, not implemented
#     "vdw-df-y"     ="sla+pw+????+vdwy"   = vdW-DF-y, reserved Thonhauser, not implemented
#     "vdw-df-z"     ="sla+pw+????+vdwz"   = vdW-DF-z, reserved Thonhauser, not implemented
#     "vdw-df-c09"   ="sla+pw+c09x+vdw1"   = vdW-DF-C09
#     "vdw-df2-c09"  ="sla+pw+c09x+vdw2"   = vdW-DF2-C09
#     "vdw-df-cx"    ="sla+pw+cx13+vdW1"   = vdW-DF-cx
#     "vdw-df-obk8"  ="sla+pw+obk8+vdw1"   = vdW-DF-obk8 (optB88-vdW)
#     "vdw-df-ob86"  ="sla+pw+ob86+vdw1"   = vdW-DF-ob86 (optB86b-vdW)
#     "vdw-df2-b86r" ="sla+pw+b86r+vdw2"   = vdW-DF2-B86R (rev-vdw-df2)
#     "rvv10" = "sla+pw+rw86+pbc+vv10"     = rVV10
#
# Note: as a rule, all keywords should be unique, and should be different
# from the short name, but there are a few exceptions.
# 
# References:
#              pz      J.P.Perdew and A.Zunger, PRB 23, 5048 (1981) 
#              vwn     S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
#              vwn1-rpa S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
#              wig     E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938) 
#              hl      L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)
#              gl      O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)
#              pw      J.P.Perdew and Y.Wang, PRB 45, 13244 (1992) 
#              obpz    G.Ortiz and P.Ballone, PRB 50, 1391 (1994) 
#              obpw    as above
#              b88     A.D.Becke, PRA 38, 3098 (1988)
#              p86     J.P.Perdew, PRB 33, 8822 (1986)
#              pw86    J.P.Perdew, PRB 33, 8800 (1986)
#              b86b    A.D.Becke, J.Chem.Phys. 85, 7184 (1986) 
#              ob86    Klimes, Bowler, Michaelides, PRB 83, 195131 (2011)
#              b86r    I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
#              pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
#              pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)
#              blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)
#              hcth    Handy et al, JCP 109, 6264 (1998)
#              olyp    Handy et al, JCP 116, 5411 (2002)
#              revPBE  Zhang and Yang, PRL 80, 890 (1998)
#              pbesol  J.P. Perdew et al., PRL 100, 136406 (2008)
#              q2d     L. Chiodo et al., PRL 108, 126402 (2012)
#              rw86    E. Amonn D. Murray et al, J. Chem. Theory comp. 5, 2754 (2009) 
#              wc      Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
#              kzk     H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)
#              pbe0    J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)
#              hse     Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
#                      Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006).
#              b3lyp   P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch
#                      J.Phys.Chem 98, 11623 (1994)
#              x3lyp   X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)
#              vdW-DF       M. Dion et al., PRL 92, 246401 (2004)
#                           T. Thonhauser et al., PRL 115, 136402 (2015)
#              vdW-DF2      Lee et al., Phys. Rev. B 82, 081101 (2010)
#              rev-vdW-DF2  I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
#              vdW-DF-cx    K. Berland and P. Hyldgaard, PRB 89, 035412 (2014)
#              vdW-DF-obk8  Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)
#              vdW-DF-ob86  Klimes et al, Phys. Rev. B, 83, 195131 (2011)
#              c09x    V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
#              tpss    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria, 
#                      PRL 91, 146401 (2003)
#              tb09    F Tran and P Blaha, Phys.Rev.Lett. 102, 226401 (2009) 
#              sogga   Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)
#              m06l    Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)
#              gau-pbe J.-W. Song, K. Yamashita, K. Hirao JCP 135, 071103 (2011)
#              rVV10   R. Sabatini et al. Phys. Rev. B 87, 041108(R) (2013)
#              ev93     Engel-Vosko, Phys. Rev. B 47, 13164 (1993)
# 
# NOTE ABOUT HSE: there are two slight deviations with respect to the HSE06 
# functional as it is in Gaussian code (that is considered as the reference
# in the chemistry community):
# - The range separation in Gaussian is precisely 0.11 bohr^-1, 
#   instead of 0.106 bohr^-1 in this implementation
# - The gradient scaling relation is a bit more complicated 
#   [ see: TM Henderson, AF Izmaylov, G Scalmani, and GE Scuseria,
#          J. Chem. Phys. 131, 044108 (2009) ]
# These two modifications accounts only for a 1e-5 Ha difference for a 
# single He atom. Info by Fabien Bruneval
