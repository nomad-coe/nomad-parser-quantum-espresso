# original developer: Henning Glawe, henning.glawe@mpsd.mpg.de
# updated by: temok.salazar@physik.hu-berlin.de

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
        try:
            xf_num_split_i = [int(x) for x in xc_functional_num.split(' ') if x.strip() != '']
        except Exception:
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
    xf_num = parse_qe_xc_num(xc_functional_num)
    LOGGER.debug('num <- input: %s <- %s, exx_fraction: %s', str(xf_num),
                 xc_functional_num, str(exact_exchange_fraction))
    # use dictionaries to ensure uniqueness:
    #   exchange/correlation functionals may be combined into _XC_, and we
    #   only want to emit such combinations once
    xc_data = {}
    #   libXC definitions of GGAs/Hybrids include the density contribution
    #   while QE allows to freely combine density, gradient etc. dependency
    #   so we need to remove at least the default settings to be compliant
    #   with libXC and NOMAD metaInfo
    xc_data_remove = {}
    # collect everything that will go to section_method
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
                xc_data, this_term, exact_exchange_fraction)
        # collect all to-be-removed terms
        for this_term in this_component.get('xc_terms_remove', []):
            apply_term_add(
                xc_data_remove, this_term, exact_exchange_fraction)
    # remove all collected terms from xc_terms_remove from xc_data
    apply_terms_remove(xc_data, xc_data_remove)
    # filter terms in xc_data:
    apply_terms_filter(xc_data)
    xc_functional = xc_functional_str(xc_data)
    # apply shortcuts/overrides for libXC compliance
    if xc_functional in LIBXC_SHORTCUT:
        xc_data = {}
        for this_term in LIBXC_SHORTCUT[xc_functional]['xc_terms']:
            apply_term_add(
                xc_data, this_term, exact_exchange_fraction)
        apply_terms_filter(xc_data)
        xc_functional = xc_functional_str(xc_data)
    xc_section_method['XC_functional'] = xc_functional
    result = []
    for k in sorted(xc_data.keys()):
        v = xc_data[k]
        result.append(v)
    return (xc_section_method, result)


def apply_term_add(xc_data, this_term, exact_exchange_fraction):
    term = copy.deepcopy(this_term)
    if 'exx_compute_weight' in term:
        term['XC_functional_weight'] = term['exx_compute_weight'](
            exact_exchange_fraction)
    if 'XC_functional_weight' not in term:
        term['XC_functional_weight'] = 1.0
    if term['XC_functional_name'] not in xc_data:
        xc_data[term['XC_functional_name']] = term
    else:
        LOGGER.info("pre-existing XC term: %s",
                    term['XC_functional_name'])
    return xc_data


def apply_terms_remove(xc_data, xc_data_remove):
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
            v.pop('exx_compute_weight', None)


def xc_functional_str(xc_data, separator='+'):
    result = ''
    for k in sorted(xc_data.keys()):
        v = xc_data[k]
        if len(result) > 0 and v.get('XC_functional_weight', 1.0) > 0:
            result += separator
        if v.get('XC_functional_weight', None) is not None:
            result += '%.3f*' % (v['XC_functional_weight'])
        result += v['XC_functional_name']
    return result











# origin: espresso-5.4.0/Modules/funct.f90
# update:
# . New exchange-correlation functionals exist in
# .     espresso-6.5.0/Modules/funct.f90
#   short comments mark the corresponding new metainfo
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
            'exx_compute_weight': lambda exx: exx,
            'XC_functional_weight': 0.25,
        }, {
            'XC_functional_name': 'LDA_X',
            'exx_compute_weight': lambda exx: (1.0 - exx),
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
            'exx_compute_weight': lambda exx: exx,
            'XC_functional_weight': 0.20,
        }, {
            'XC_functional_name': 'LDA_X',
            'exx_compute_weight': lambda exx: (1.0 - exx),
            'XC_functional_weight': 0.8,
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       "b3lp",
            'x_qe_xc_iexch_comment':    "B3LYP(Slater*0.80+HF*0.20)",
            'x_qe_xc_iexch':      7,
        },
    },
    # LDA_X_KZK is not part of libXC. Look up it at
    # 'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional'
    {
        'xc_terms': [{
            'XC_functional_name': "LDA_X_KZK",
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
            'exx_compute_weight': lambda exx: exx,
            'XC_functional_weight': 0.218,
        }, {
            'XC_functional_name': 'LDA_X',
            'exx_compute_weight': lambda exx: (1.0 - exx),
            'XC_functional_weight': 0.782,
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name':       "x3lp",
            'x_qe_xc_iexch_comment':    "X3LYP(Slater*0.782+HF*0.218)",
            'x_qe_xc_iexch':      9,
        },
    },
    # update for espresso-6.5.0: KLI
        {
        'xc_terms': [{
            'XC_functional_name': "LDA_X_KLI",
        }],
        'xc_section_method': {
            'x_qe_xc_iexch_name': "kli",
            'x_qe_xc_iexch_comment': "KLI aproximation for exx",
            'x_qe_xc_iexch': 10,
        },
    },
]

# Correlation functionals UNchanged between espresso v5.4 & v6.5
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
            'XC_functional_name': "LDA_C_KZK",
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
            'XC_functional_name': 'LDA_C_VWN_RPA',
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
            'XC_functional_name': 'LDA_C_VWN_RPA',
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

# New "exchange_gradient_correction" functionals for q-espresso (qe) v6.5
#    igcx=[1..28] unchanged between qe-v5.4 & v6.5
# New additions:  igcx=[29..42]

EXCHANGE_GRADIENT_CORRECTION = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B88",
        }],
        'xc_terms_remove': [{
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
        'xc_terms_remove': [{
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
        'xc_terms_remove': [{
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
        'xc_terms_remove': [{
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
        'xc_terms_remove': [{
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
            'exx_compute_weight': lambda exx: (1.0 - exx),
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
            'XC_functional_weight': 0.75,
            'exx_compute_weight': lambda exx: (1.0 - exx),
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
            'exx_compute_weight': lambda exx: 0.72 if abs(exx) > 0.01 else 1.0
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
            'XC_functional_weight': 0.8,
            'exx_compute_weight': lambda exx: (1.0 - exx),
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
            'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > 0.01) else 0.0
        }, {
            'XC_functional_name': "GGA_X_PBE",
            'exx_compute_weight': lambda exx: 0.0 if (abs(exx) > 0.01) else 1.0
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }, {
            'XC_functional_name': 'GGA_C_PBE',
            'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > 0.01) else 0.0
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "q2dx",
            'x_qe_xc_igcx_comment':    "Q2D exchange grad corr",
            'x_qe_xc_igcx':      19,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "HYB_GGA_XC_GAU_PBE",
            'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > 0.01) else 0.0
        }, {
            'XC_functional_name': "GGA_X_PBE",
            'exx_compute_weight': lambda exx: 0.0 if (abs(exx) > 0.01) else 1.0
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }, {
            'XC_functional_name': 'GGA_C_PBE',
            'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > 0.01) else 0.0
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':       "obk8",
            'x_qe_xc_igcx_comment':    "optB88  exchange",
            'x_qe_xc_igcx':      23,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_OPTB86_VDW",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
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
            'exx_compute_weight':
                lambda exx: 0.542 if (abs(exx) > 0.01) else 1.0
        }, {
            'XC_functional_name': "GGA_X_PW91",
            'XC_functional_weight': 0.167,
            'exx_compute_weight':
                lambda exx: 0.167 if (abs(exx) > 0.01) else 0.0
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
            'exx_compute_weight':
                lambda exx: 0.709 if (abs(exx) > 0.01) else 1.0
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':    "x3lp",
            'x_qe_xc_igcx_comment': "X3LYP (Becke88*0.542 "
                                    " + Perdew-Wang91*0.167)",
            'x_qe_xc_igcx':         28,
        },
    },
# New additions for qe-v6.5.0:  igcx=[29..42]
# - - - - - -
# igcx: 29. The ingredient 'vdW-DF-cx' is documented in the nomad-meta-info, where it
# has the name 'vdw_c_df_cx'
# 'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/-/wikis/metainfo/XC-functional':
# 'vdW-DF-cx' implies igcx=27 => 'cx13', hence LDA_X is implicit. Full weight.
    {
        'xc_terms': [{
            'XC_functional_name': "vdw_c_df_cx",
        }, {
        'XC_functional_name': 'HF_X',
        'exx_compute_weight': lambda exx: exx,
        'XC_functional_weight': 0.25,
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':    "cx0",
            'x_qe_xc_igcx_comment': "vdW-DF-cx+HF/4 (cx13-0)",
            'x_qe_xc_igcx':         29,
        },
    },
# - - - - - -
# igcx:30. Needs full LDA_X removal, due to 'GGA_X_RPW86' (see igcx:27)
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RPW86",
        },{
            'XC_functional_name': 'HF_X',
            'exx_compute_weight': lambda exx: exx,
            'XC_functional_weight': 0.25,
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "r860",
            'x_qe_xc_igcx_comment': "rPW86+HF/4 (rw86-0); (for DF0)",
            'x_qe_xc_igcx': 30,
        },
    },
# - - - - - -
# igcx:31. Similar comments as in 'igcx:29'
    {
        'xc_terms': [{
            'XC_functional_name': "vdw_c_df_cx",
        }, {
            'XC_functional_name': 'HF_X',
            'exx_compute_weight': lambda exx: exx,
            'XC_functional_weight': 0.20,
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name':    "cx0p",
            'x_qe_xc_igcx_comment': "vdW-DF-cx+HF/5 (cx13-0p)",
            'x_qe_xc_igcx':         31,
        },
    },
# - - - - - -
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RESERVED",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "ahcx",
            'x_qe_xc_igcx_comment': "vdW-DF-cx based; not yet in use (reserved PH)",
            'x_qe_xc_igcx': 32,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RESERVED",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "ahf2",
            'x_qe_xc_igcx_comment': "vdW-DF2 based; not yet in use (reserved PH)",
            'x_qe_xc_igcx': 33,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RESERVED",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "ahpb",
            'x_qe_xc_igcx_comment': "PBE based; not yet in use (reserved PH)",
            'x_qe_xc_igcx': 34,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RESERVED",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "ahps",
            'x_qe_xc_igcx_comment': "PBE-sol based; not in use (reserved PH)",
            'x_qe_xc_igcx': 35,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RESERVED",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "cx14",
            'x_qe_xc_igcx_comment': "Exporations (typo?: explorations), (reserved PH)",
            'x_qe_xc_igcx': 36,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RESERVED",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "cx15",
            'x_qe_xc_igcx_comment': "Exporations (typo? explorations?)(reserved PH)",
            'x_qe_xc_igcx': 37,
        },
    },
#
# igcx': 38. Ingredients:
#   'b86r' -> 'igcx:26' -> 'GGA_X_B86_R'
#   'vdW-DF2' -> 'vdw_c_df2' . See nomad's gitlab:
#   'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/-/wikis/metainfo/XC-functional':

    {
     'xc_terms': [{
            'XC_functional_name': "vdw_c_df2",
        }, {
            'XC_functional_name': "GGA_X_B86_R",
        }, {
            'XC_functional_name': 'HF_X',
            'exx_compute_weight': lambda exx: exx,
            'XC_functional_weight': 0.25,
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "br0",
            'x_qe_xc_igcx_comment': "vdW-DF2-b86r+HF/4 (b86r-0)",
            'x_qe_xc_igcx': 38,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_RESERVED",
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "cx16",
            'x_qe_xc_igcx_comment': "Exporations (typo?, explorations?)(reserved PH)",
            'x_qe_xc_igcx': 39,
        },
    },
#- - - -
    {
        'xc_terms': [{
            'XC_functional_name': "vdw_c_df1",
        }, {
            'XC_functional_name': "GGA_X_C09X",
        }, {
            'XC_functional_name': 'HF_X',
            'exx_compute_weight': lambda exx: exx,
            'XC_functional_weight': 0.25,
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "c090",
            'x_qe_xc_igcx_comment': "vdW-DF-c09+HF/4 (c09-0)",
            'x_qe_xc_igcx': 40,
        },
    },
# - - - - - - -
# 'igcx:41' Note: 'B86b' is defined in 'igcx:22'
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B86_MGC",
        }],
            'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
            'XC_functional_weight': 0.75,
            'exx_compute_weight': lambda exx: (1.0 - exx),
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "b86x",
            'x_qe_xc_igcx_comment': "B86b exchange * 0.75",
            'x_qe_xc_igcx': 41,
        },
    },
# - - - - - - -
# 'B88' is defined in 'igcx:1'
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_X_B88",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
            'XC_functional_weight': 0.50,
            'exx_compute_weight': lambda exx: (1.0 - exx),
        }],
        'xc_section_method': {
            'x_qe_xc_igcx_name': "b88x",
            'x_qe_xc_igcx_comment': "B88 exchange * 0.50",
            'x_qe_xc_igcx': 42,
        },
    },
]

# UNchanged between espresso v5.4 & v6.5
CORRELATION_GRADIENT_CORRECTION = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_P86",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_LYP',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_LYP',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
            'x_qe_xc_igcc_name':    "m6lc",
            'x_qe_xc_igcc_comment': "M06L Meta-GGA (Espresso-version < 5.1)",
            'x_qe_xc_igcc':         11,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "GGA_C_Q2D",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_LYP',
            'XC_functional_weight': 0.871,
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':       "x3lp",
            'x_qe_xc_igcc_comment':    "X3LYP (Lee-Yang-Parr*0.871)",
            'x_qe_xc_igcc':      13,
        },
    },
    {
    #  igcc=14 is not defined in NEITHER of v5.1, v6.1, v6.4's Modules/funct.f90
    # "BEEF-vdW, a GGA with vdW-DF2 type nonlocal correlation"
        'xc_terms': [{
            'XC_functional_name':  "GGA_C_BEEF-vdW",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_C_PW',
        }],
        'xc_section_method': {
            'x_qe_xc_igcc_name':  "BEEF-vdW",
            'x_qe_xc_igcc_comment':  "libbeef V0.1.1 library",
            'x_qe_xc_igcc':  14,
        },
    },
]

# New additions for espresso-6.5.0:  imeta=[4, 5, 6]
META_GGA = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_TPSS",
        }, {
            'XC_functional_name': "MGGA_C_TPSS",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }, {
            'XC_functional_name': 'LDA_C_PW',
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
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }, {
            'XC_functional_name': 'LDA_C_PW',
        }],
        'xc_section_method': {
            'x_qe_xc_imeta_name':       "tb09",
            'x_qe_xc_imeta_comment':    "TB09 Meta-GGA",
            'x_qe_xc_imeta':      3,
        },
    },
    # imeta = [4,5,6] are new espresso-6.5.0/Modules/funct.f90
    {
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_TPSS",
        }, {
            'XC_functional_name': "MGGA_C_TPSS",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }, {
            'XC_functional_name': 'LDA_C_PW',
        }],
        'xc_section_method': {
            'x_qe_xc_imeta_name': "+meta",
            'x_qe_xc_imeta_comment': "activate MGGA even without MGGA-XC",
            'x_qe_xc_imeta': 4,
        },
    },
    # - - - - - - - -
    {
        'xc_terms': [{
            'XC_functional_name': "MGGA_X_SCAN",
        }, {
            'XC_functional_name': "MGGA_C_SCAN",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }, {
            'XC_functional_name': 'LDA_C_PW',
        }],
        'xc_section_method': {
            'x_qe_xc_imeta_name': "scan",
            'x_qe_xc_imeta_comment': "SCAN Meta-GGA ",
            'x_qe_xc_imeta': 5,
        },
    },
# - - - - - - - -
    {
        'xc_terms': [{
            'XC_functional_name': "HYB_MGGA_X_SCAN0",
        }, {
            'XC_functional_name': "MGGA_C_SCAN",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': 'LDA_X',
        }, {
            'XC_functional_name': 'LDA_C_PW',
        }],
        'xc_section_method': {
            'x_qe_xc_imeta_name': "sca0",
            'x_qe_xc_imeta_comment': "SCAN0  Meta-GGA",
            'x_qe_xc_imeta': 6,
        },
    },
]

VAN_DER_WAALS = [
    None,
    {
        'xc_terms': [{
            'XC_functional_name': "VDW_C_DF1",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': "LDA_C_PW",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vdw1",
            'x_qe_xc_inlc_comment':    "vdW-DF1",
            'x_qe_xc_inlc':      1,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "VDW_C_DF2",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': "LDA_C_PW",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vdw2",
            'x_qe_xc_inlc_comment':    "vdW-DF2",
            'x_qe_xc_inlc':      2,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "VDW_C_RVV10",
        }],
        'xc_terms_remove': [{
            'XC_functional_name': "GGA_C_PBE",
        }],
        'xc_section_method': {
            'x_qe_xc_inlc_name':       "vv10",
            'x_qe_xc_inlc_comment':    "rVV10",
            'x_qe_xc_inlc':      3,
        },
    },
    {
        'xc_terms': [{
            'XC_functional_name': "VDW_DFX_x_qe",
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
            'XC_functional_name': "VDW_DFY_x_qe",
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
            'XC_functional_name': "VDW_DFZ_x_qe",
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

LIBXC_SHORTCUT = {
    "0.810*GGA_C_LYP+0.720*GGA_X_B88+0.200*HF_X+0.190*LDA_C_VWN": {
        'xc_terms': [{
            'XC_functional_name': "HYB_GGA_XC_B3LYP",
        }]
    },
    "GGA_C_PBE+0.750*GGA_X_PBE+0.250*HF_X": {
        'xc_terms': [{
            'XC_functional_name': "HYB_GGA_XC_PBEH",
        }]
    },
}
