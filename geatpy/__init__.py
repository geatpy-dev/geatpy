# -*- coding: utf-8 -*-
"""
import all libs of geatpy

"""

import sys

__author__ = "Geatpy Team"
__version__ = "2.6.0"

# import the core
lib_path = __file__[:-11] + 'core/'
if lib_path not in sys.path:
    sys.path.append(lib_path)
from awGA import awGA
from boundfix import boundfix
from bs2int import bs2int
from bs2real import bs2real
from bs2ri import bs2ri
from crowdis import crowdis
from crtbp import crtbp
from crtfld import crtfld
from crtgp import crtgp
from crtidp import crtidp
from crtip import crtip
from crtpc import crtpc
from crtpp import crtpp
from crtri import crtri
from crtrp import crtrp
from crtup import crtup
from dup import dup
from ecs import ecs
from etour import etour
from indexing import indexing
import indicator
from mergecv import mergecv
from migrate import migrate
from moeaplot import moeaplot
from mselecting import mselecting
from mutate import mutate
from mutbga import mutbga
from mutbin import mutbin
from mutde import mutde
from mutgau import mutgau
from mutinv import mutinv
from mutmove import mutmove
from mutpolyn import mutpolyn
from mutpp import mutpp
from mutswap import mutswap
from mutuni import mutuni
from ndsortDED import ndsortDED
from ndsortESS import ndsortESS
from ndsortTNS import ndsortTNS
from otos import otos
from pbi import pbi
from powing import powing
from ranking import ranking
from rcs import rcs
from recdis import recdis
from recint import recint
from reclin import reclin
from recndx import recndx
from recombin import recombin
from recsbx import recsbx
from refgselect import refgselect
from refselect import refselect
from ri2bs import ri2bs
from rps import rps
from rwGA import rwGA
from rws import rws
from scaling import scaling
from selecting import selecting
from soeaplot import soeaplot
from sus import sus
from tcheby import tcheby
from tour import tour
from trcplot import trcplot
from urs import urs
from varplot import varplot
from xovbd import xovbd
from xovdp import xovdp
from xovexp import xovexp
from xovox import xovox
from xovpmx import xovpmx
from xovsec import xovsec
from xovsh import xovsh
from xovsp import xovsp
from xovud import xovud

# ========================import the EA framework=======================
# import Classes and mathods
lib_path = __file__[:-11]
if lib_path not in sys.path:
    sys.path.append(lib_path)
from Algorithm import Algorithm
from Algorithm import MoeaAlgorithm
from Algorithm import SoeaAlgorithm
from Population import Population
from Problem import Problem
from PsyPopulation import PsyPopulation

# import operators
from operators.migration.Migrate import Migrate
from operators.mutation.Mutation import Mutation
from operators.mutation.Mutbga import Mutbga
from operators.mutation.Mutbin import Mutbin
from operators.mutation.Mutde import Mutde
from operators.mutation.Mutgau import Mutgau
from operators.mutation.Mutinv import Mutinv
from operators.mutation.Mutmove import Mutmove
from operators.mutation.Mutpolyn import Mutpolyn
from operators.mutation.Mutpp import Mutpp
from operators.mutation.Mutswap import Mutswap
from operators.mutation.Mutuni import Mutuni
from operators.recombination.Recombination import Recombination
from operators.recombination.Recdis import Recdis
from operators.recombination.Recint import Recint
from operators.recombination.Reclin import Reclin
from operators.recombination.Recndx import Recndx
from operators.recombination.Recsbx import Recsbx
from operators.recombination.Xovbd import Xovbd
from operators.recombination.Xovdp import Xovdp
from operators.recombination.Xovexp import Xovexp
from operators.recombination.Xovox import Xovox
from operators.recombination.Xovpmx import Xovpmx
from operators.recombination.Xovsec import Xovsec
from operators.recombination.Xovsh import Xovsh
from operators.recombination.Xovsp import Xovsp
from operators.recombination.Xovud import Xovud

# import soea templates
from templates.soeas.DE.DE_best_1_bin.soea_DE_best_1_bin_templet import soea_DE_best_1_bin_templet
from templates.soeas.DE.DE_best_1_L.soea_DE_best_1_L_templet import soea_DE_best_1_L_templet
from templates.soeas.DE.DE_currentToBest_1_bin.soea_DE_currentToBest_1_bin_templet import \
    soea_DE_currentToBest_1_bin_templet
from templates.soeas.DE.DE_currentToBest_1_L.soea_DE_currentToBest_1_L_templet import soea_DE_currentToBest_1_L_templet
from templates.soeas.DE.DE_currentToRand_1.soea_DE_currentToRand_1_templet import soea_DE_currentToRand_1_templet
from templates.soeas.DE.DE_rand_1_bin.soea_DE_rand_1_bin_templet import soea_DE_rand_1_bin_templet
from templates.soeas.DE.DE_rand_1_L.soea_DE_rand_1_L_templet import soea_DE_rand_1_L_templet
from templates.soeas.DE.DE_targetToBest_1_bin.soea_DE_targetToBest_1_bin_templet import \
    soea_DE_targetToBest_1_bin_templet
from templates.soeas.DE.DE_targetToBest_1_L.soea_DE_targetToBest_1_L_templet import soea_DE_targetToBest_1_L_templet
from templates.soeas.ES.ES_1_plus_1_templet.soea_ES_1_plus_1_templet import soea_ES_1_plus_1_templet
from templates.soeas.ES.ES_miu_plus_lambda_templet.soea_ES_miu_plus_lambda_templet import \
    soea_ES_miu_plus_lambda_templet
from templates.soeas.GA.EGA.soea_psy_EGA_templet import soea_psy_EGA_templet
from templates.soeas.GA.EGA.soea_EGA_templet import soea_EGA_templet
from templates.soeas.GA.SEGA.soea_multi_SEGA_templet import soea_multi_SEGA_templet
from templates.soeas.GA.SEGA.soea_psy_SEGA_templet import soea_psy_SEGA_templet
from templates.soeas.GA.SEGA.soea_SEGA_templet import soea_SEGA_templet
from templates.soeas.GA.SGA.soea_GGAP_SGA_templet import soea_GGAP_SGA_templet
from templates.soeas.GA.SGA.soea_psy_GGAP_SGA_templet import soea_psy_GGAP_SGA_templet
from templates.soeas.GA.SGA.soea_psy_SGA_templet import soea_psy_SGA_templet
from templates.soeas.GA.SGA.soea_SGA_templet import soea_SGA_templet
from templates.soeas.GA.steadyGA.soea_psy_steadyGA_templet import soea_psy_steadyGA_templet
from templates.soeas.GA.steadyGA.soea_steadyGA_templet import soea_steadyGA_templet
from templates.soeas.GA.studGA.soea_psy_studGA_templet import soea_psy_studGA_templet
from templates.soeas.GA.studGA.soea_studGA_templet import soea_studGA_templet

# import moea templates
from templates.moeas.awGA.moea_awGA_templet import moea_awGA_templet
from templates.moeas.awGA.moea_psy_awGA_templet import moea_psy_awGA_templet
from templates.moeas.moead.moea_MOEAD_archive_templet import moea_MOEAD_archive_templet
from templates.moeas.moead.moea_MOEAD_templet import moea_MOEAD_templet
from templates.moeas.moead.moea_MOEAD_DE_templet import moea_MOEAD_DE_templet
from templates.moeas.nsga2.moea_NSGA2_archive_templet import moea_NSGA2_archive_templet
from templates.moeas.nsga2.moea_NSGA2_DE_templet import moea_NSGA2_DE_templet
from templates.moeas.nsga2.moea_NSGA2_templet import moea_NSGA2_templet
from templates.moeas.nsga2.moea_psy_NSGA2_archive_templet import moea_psy_NSGA2_archive_templet
from templates.moeas.nsga2.moea_psy_NSGA2_templet import moea_psy_NSGA2_templet
from templates.moeas.nsga3.moea_NSGA3_DE_templet import moea_NSGA3_DE_templet
from templates.moeas.nsga3.moea_NSGA3_templet import moea_NSGA3_templet
from templates.moeas.nsga3.moea_psy_NSGA3_templet import moea_psy_NSGA3_templet
from templates.moeas.pps.moea_PPS_MOEAD_DE_archive_templet import moea_PPS_MOEAD_DE_archive_templet
from templates.moeas.rvea.moea_psy_RVEA_RES_templet import moea_psy_RVEA_RES_templet
from templates.moeas.rvea.moea_psy_RVEA_templet import moea_psy_RVEA_templet
from templates.moeas.rvea.moea_RVEA_templet import moea_RVEA_templet
from templates.moeas.rvea.moea_RVEA_RES_templet import moea_RVEA_RES_templet
