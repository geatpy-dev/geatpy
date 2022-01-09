# -*- coding: utf-8 -*-
"""
import all libs of geatpy

"""

__author__ = "Geatpy Team"
__version__ = "2.7.0"

# import the core
from geatpy.core.awGA import awGA
from geatpy.core.boundfix import boundfix
from geatpy.core.bs2int import bs2int
from geatpy.core.bs2real import bs2real
from geatpy.core.bs2ri import bs2ri
from geatpy.core.cdist import cdist
from geatpy.core.crowdis import crowdis
from geatpy.core.crtbp import crtbp
from geatpy.core.crtfld import crtfld
from geatpy.core.crtgp import crtgp
from geatpy.core.crtidp import crtidp
from geatpy.core.crtip import crtip
from geatpy.core.crtpc import crtpc
from geatpy.core.crtpp import crtpp
from geatpy.core.crtri import crtri
from geatpy.core.crtrp import crtrp
from geatpy.core.crtup import crtup
from geatpy.core.dup import dup
from geatpy.core.ecs import ecs
from geatpy.core.etour import etour
from geatpy.core.indexing import indexing
import geatpy.core.indicator
from geatpy.core.mergecv import mergecv
from geatpy.core.migrate import migrate
from geatpy.core.moeaplot import moeaplot
from geatpy.core.mselecting import mselecting
from geatpy.core.mutate import mutate
from geatpy.core.mutbga import mutbga
from geatpy.core.mutbin import mutbin
from geatpy.core.mutde import mutde
from geatpy.core.mutgau import mutgau
from geatpy.core.mutinv import mutinv
from geatpy.core.mutmove import mutmove
from geatpy.core.mutpolyn import mutpolyn
from geatpy.core.mutpp import mutpp
from geatpy.core.mutswap import mutswap
from geatpy.core.mutuni import mutuni
from geatpy.core.ndsortDED import ndsortDED
from geatpy.core.ndsortESS import ndsortESS
from geatpy.core.ndsortTNS import ndsortTNS
from geatpy.core.otos import otos
from geatpy.core.pbi import pbi
from geatpy.core.powing import powing
from geatpy.core.ranking import ranking
from geatpy.core.rcs import rcs
from geatpy.core.recdis import recdis
from geatpy.core.recint import recint
from geatpy.core.reclin import reclin
from geatpy.core.recndx import recndx
from geatpy.core.recombin import recombin
from geatpy.core.recsbx import recsbx
from geatpy.core.refgselect import refgselect
from geatpy.core.refselect import refselect
from geatpy.core.ri2bs import ri2bs
from geatpy.core.rps import rps
from geatpy.core.rwGA import rwGA
from geatpy.core.rws import rws
from geatpy.core.scaling import scaling
from geatpy.core.selecting import selecting
from geatpy.core.soeaplot import soeaplot
from geatpy.core.sus import sus
from geatpy.core.tcheby import tcheby
from geatpy.core.tour import tour
from geatpy.core.trcplot import trcplot
from geatpy.core.urs import urs
from geatpy.core.varplot import varplot
from geatpy.core.xovbd import xovbd
from geatpy.core.xovdp import xovdp
from geatpy.core.xovexp import xovexp
from geatpy.core.xovox import xovox
from geatpy.core.xovpmx import xovpmx
from geatpy.core.xovsec import xovsec
from geatpy.core.xovsh import xovsh
from geatpy.core.xovsp import xovsp
from geatpy.core.xovud import xovud

# import classes and mathods
from geatpy.Algorithm import Algorithm
from geatpy.Algorithm import MoeaAlgorithm
from geatpy.Algorithm import SoeaAlgorithm
from geatpy.optimize import optimize
from geatpy.Population import Population
from geatpy.Problem import Problem
from geatpy.PsyPopulation import PsyPopulation
from geatpy.core import indicator

# import visualization
from geatpy.visualization.PointScatter import PointScatter
from geatpy.visualization.ParCoordPlotter import ParCoordPlotter

# import benchmarks
import geatpy.benchmarks

# import operators
from geatpy.operators.migration.Migrate import Migrate
from geatpy.operators.mutation.Mutation import Mutation
from geatpy.operators.mutation.Mutbga import Mutbga
from geatpy.operators.mutation.Mutbin import Mutbin
from geatpy.operators.mutation.Mutde import Mutde
from geatpy.operators.mutation.Mutgau import Mutgau
from geatpy.operators.mutation.Mutinv import Mutinv
from geatpy.operators.mutation.Mutmove import Mutmove
from geatpy.operators.mutation.Mutpolyn import Mutpolyn
from geatpy.operators.mutation.Mutpp import Mutpp
from geatpy.operators.mutation.Mutswap import Mutswap
from geatpy.operators.mutation.Mutuni import Mutuni
from geatpy.operators.recombination.Recdis import Recdis
from geatpy.operators.recombination.Recint import Recint
from geatpy.operators.recombination.Reclin import Reclin
from geatpy.operators.recombination.Recndx import Recndx
from geatpy.operators.recombination.Recombination import Recombination
from geatpy.operators.recombination.Recsbx import Recsbx
from geatpy.operators.recombination.Xovbd import Xovbd
from geatpy.operators.recombination.Xovdp import Xovdp
from geatpy.operators.recombination.Xovexp import Xovexp
from geatpy.operators.recombination.Xovox import Xovox
from geatpy.operators.recombination.Xovpmx import Xovpmx
from geatpy.operators.recombination.Xovsec import Xovsec
from geatpy.operators.recombination.Xovsh import Xovsh
from geatpy.operators.recombination.Xovsp import Xovsp
from geatpy.operators.recombination.Xovud import Xovud
# import moea algorithms
from geatpy.algorithms.moeas.awGA.moea_awGA_templet import moea_awGA_templet
from geatpy.algorithms.moeas.awGA.moea_psy_awGA_templet import moea_psy_awGA_templet
from geatpy.algorithms.moeas.moead.moea_MOEAD_DE_templet import moea_MOEAD_DE_templet
from geatpy.algorithms.moeas.moead.moea_MOEAD_archive_templet import moea_MOEAD_archive_templet
from geatpy.algorithms.moeas.moead.moea_MOEAD_templet import moea_MOEAD_templet
from geatpy.algorithms.moeas.nsga2.moea_NSGA2_DE_templet import moea_NSGA2_DE_templet
from geatpy.algorithms.moeas.nsga2.moea_NSGA2_archive_templet import moea_NSGA2_archive_templet
from geatpy.algorithms.moeas.nsga2.moea_NSGA2_templet import moea_NSGA2_templet
from geatpy.algorithms.moeas.nsga2.moea_psy_NSGA2_archive_templet import moea_psy_NSGA2_archive_templet
from geatpy.algorithms.moeas.nsga2.moea_psy_NSGA2_templet import moea_psy_NSGA2_templet
from geatpy.algorithms.moeas.nsga3.moea_NSGA3_DE_templet import moea_NSGA3_DE_templet
from geatpy.algorithms.moeas.nsga3.moea_NSGA3_templet import moea_NSGA3_templet
from geatpy.algorithms.moeas.nsga3.moea_psy_NSGA3_templet import moea_psy_NSGA3_templet
from geatpy.algorithms.moeas.pps.moea_PPS_MOEAD_DE_archive_templet import moea_PPS_MOEAD_DE_archive_templet
from geatpy.algorithms.moeas.rvea.moea_RVEA_RES_templet import moea_RVEA_RES_templet
from geatpy.algorithms.moeas.rvea.moea_RVEA_templet import moea_RVEA_templet
from geatpy.algorithms.moeas.rvea.moea_psy_RVEA_RES_templet import moea_psy_RVEA_RES_templet
from geatpy.algorithms.moeas.rvea.moea_psy_RVEA_templet import moea_psy_RVEA_templet
from geatpy.algorithms.soeas.DE.DE_best_1_L.soea_DE_best_1_L_templet import soea_DE_best_1_L_templet
# import soea algorithms
from geatpy.algorithms.soeas.DE.DE_best_1_bin.soea_DE_best_1_bin_templet import soea_DE_best_1_bin_templet
from geatpy.algorithms.soeas.DE.DE_currentToBest_1_L.soea_DE_currentToBest_1_L_templet import \
    soea_DE_currentToBest_1_L_templet
from geatpy.algorithms.soeas.DE.DE_currentToBest_1_bin.soea_DE_currentToBest_1_bin_templet import \
    soea_DE_currentToBest_1_bin_templet
from geatpy.algorithms.soeas.DE.DE_currentToRand_1.soea_DE_currentToRand_1_templet import soea_DE_currentToRand_1_templet
from geatpy.algorithms.soeas.DE.DE_rand_1_L.soea_DE_rand_1_L_templet import soea_DE_rand_1_L_templet
from geatpy.algorithms.soeas.DE.DE_rand_1_bin.soea_DE_rand_1_bin_templet import soea_DE_rand_1_bin_templet
from geatpy.algorithms.soeas.DE.DE_targetToBest_1_L.soea_DE_targetToBest_1_L_templet import \
    soea_DE_targetToBest_1_L_templet
from geatpy.algorithms.soeas.DE.DE_targetToBest_1_bin.soea_DE_targetToBest_1_bin_templet import \
    soea_DE_targetToBest_1_bin_templet
from geatpy.algorithms.soeas.ES.ES_1_plus_1.soea_ES_1_plus_1_templet import soea_ES_1_plus_1_templet
from geatpy.algorithms.soeas.ES.ES_miu_plus_lambda.soea_ES_miu_plus_lambda_templet import \
    soea_ES_miu_plus_lambda_templet
from geatpy.algorithms.soeas.GA.EGA.soea_EGA_templet import soea_EGA_templet
from geatpy.algorithms.soeas.GA.EGA.soea_psy_EGA_templet import soea_psy_EGA_templet
from geatpy.algorithms.soeas.GA.SEGA.soea_SEGA_templet import soea_SEGA_templet
from geatpy.algorithms.soeas.GA.SEGA.soea_multi_SEGA_templet import soea_multi_SEGA_templet
from geatpy.algorithms.soeas.GA.SEGA.soea_psy_SEGA_templet import soea_psy_SEGA_templet
from geatpy.algorithms.soeas.GA.SGA.soea_GGAP_SGA_templet import soea_GGAP_SGA_templet
from geatpy.algorithms.soeas.GA.SGA.soea_SGA_templet import soea_SGA_templet
from geatpy.algorithms.soeas.GA.SGA.soea_psy_GGAP_SGA_templet import soea_psy_GGAP_SGA_templet
from geatpy.algorithms.soeas.GA.SGA.soea_psy_SGA_templet import soea_psy_SGA_templet
from geatpy.algorithms.soeas.GA.steadyGA.soea_psy_steadyGA_templet import soea_psy_steadyGA_templet
from geatpy.algorithms.soeas.GA.steadyGA.soea_steadyGA_templet import soea_steadyGA_templet
from geatpy.algorithms.soeas.GA.studGA.soea_psy_studGA_templet import soea_psy_studGA_templet
from geatpy.algorithms.soeas.GA.studGA.soea_studGA_templet import soea_studGA_templet
