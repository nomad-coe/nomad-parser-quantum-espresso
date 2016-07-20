import logging
import re


LOGGER = logging.getLogger(__name__)


# origin: espresso-5.4.0/Modules/funct.f90
# Any nonconflicting combination of the following keywords is acceptable:
# 
# Exchange:    "nox"    none                           iexch=0
#              "sla"    Slater (alpha=2/3)             iexch=1 (default)
#              "sl1"    Slater (alpha=1.0)             iexch=2
#              "rxc"    Relativistic Slater            iexch=3
#              "oep"    Optimized Effective Potential  iexch=4
#              "hf"     Hartree-Fock                   iexch=5
#              "pb0x"   PBE0 (Slater*0.75+HF*0.25)     iexch=6
#              "b3lp"   B3LYP(Slater*0.80+HF*0.20)     iexch=7
#              "kzk"    Finite-size corrections        iexch=8
#              "x3lp"   X3LYP(Slater*0.782+HF*0.218)   iexch=9
# 
# Correlation: "noc"    none                           icorr=0
#              "pz"     Perdew-Zunger                  icorr=1 (default)
#              "vwn"    Vosko-Wilk-Nusair              icorr=2
#              "lyp"    Lee-Yang-Parr                  icorr=3
#              "pw"     Perdew-Wang                    icorr=4
#              "wig"    Wigner                         icorr=5
#              "hl"     Hedin-Lunqvist                 icorr=6
#              "obz"    Ortiz-Ballone form for PZ      icorr=7
#              "obw"    Ortiz-Ballone form for PW      icorr=8
#              "gl"     Gunnarson-Lunqvist             icorr=9
#              "kzk"    Finite-size corrections        icorr=10
#              "vwn-rpa" Vosko-Wilk-Nusair, alt param  icorr=11
#              "b3lp"   B3LYP (0.19*vwn+0.81*lyp)      icorr=12
#              "b3lpv1r"  B3LYP-VWN-1-RPA 
#                         (0.19*vwn_rpa+0.81*lyp)      icorr=13
#              "x3lp"   X3LYP (0.129*vwn_rpa+0.871*lyp)icorr=14
# 
# Gradient Correction on Exchange:
#              "nogx"   none                           igcx =0 (default)
#              "b88"    Becke88 (beta=0.0042)          igcx =1
#              "ggx"    Perdew-Wang 91                 igcx =2
#              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
#              "rpb"    revised PBE by Zhang-Yang      igcx =4
#              "hcth"   Cambridge exch, Handy et al    igcx =5
#              "optx"   Handy's exchange functional    igcx =6
#              "pb0x"   PBE0 (PBE exchange*0.75)       igcx =8
#              "b3lp"   B3LYP (Becke88*0.72)           igcx =9
#              "psx"    PBEsol exchange                igcx =10
#              "wcx"    Wu-Cohen                       igcx =11
#              "hse"    HSE screened exchange          igcx =12
#              "rw86"   revised PW86                   igcx =13
#              "pbe"    same as PBX, back-comp.        igcx =14
#              "c09x"   Cooper 09                      igcx =16
#              "sox"    sogga                          igcx =17
#              "q2dx"   Q2D exchange grad corr         igcx =19
#              "gaup"   Gau-PBE hybrid exchange        igcx =20
#              "pw86"   Perdew-Wang (1986) exchange    igcx =21
#              "b86b"   Becke (1986) exchange          igcx =22
#              "obk8"   optB88  exchange               igcx =23
#              "ob86"   optB86b exchange               igcx =24
#              "evx"    Engel-Vosko exchange           igcx =25
#              "b86r"   revised Becke (b86b)           igcx =26
#              "cx13"   consistent exchange            igcx =27
#              "x3lp"   X3LYP (Becke88*0.542 +
#                              Perdew-Wang91*0.167)    igcx =28
# 
# Gradient Correction on Correlation:
#              "nogc"   none                           igcc =0 (default)
#              "p86"    Perdew86                       igcc =1
#              "ggc"    Perdew-Wang 91 corr.           igcc =2
#              "blyp"   Lee-Yang-Parr                  igcc =3
#              "pbc"    Perdew-Burke-Ernzenhof corr    igcc =4
#              "hcth"   Cambridge corr, Handy et al    igcc =5
#              "b3lp"   B3LYP (Lee-Yang-Parr*0.81)     igcc =7
#              "psc"    PBEsol corr                    igcc =8
#              "pbe"    same as PBX, back-comp.        igcc =9
#              "q2dc"   Q2D correlation grad corr      igcc =12
#              "x3lp"   X3LYP (Lee-Yang-Parr*0.871)    igcc =13
# 
# Meta-GGA functionals
#              "tpss"   TPSS Meta-GGA                  imeta=1
#              "m6lx"   M06L Meta-GGA                  imeta=2
#              "tb09"   TB09 Meta-GGA                  imeta=3
# 
# Van der Waals functionals (nonlocal term only)
#              "nonlc"  none                           inlc =0 (default)
#              "vdw1"   vdW-DF1                        inlc =1
#              "vdw2"   vdW-DF2                        inlc =2
#              "vv10"   rVV10                          inlc =3
#              "vdwx"   vdW-DF-x                       inlc =4, reserved Thonhauser, not implemented
#              "vdwy"   vdW-DF-y                       inlc =5, reserved Thonhauser, not implemented
#              "vdwz"   vdW-DF-z                       inlc =6, reserved Thonhauser, not implemented
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

