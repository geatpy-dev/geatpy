""" It's a test of importing geatpy"""

import sys
import platform

lib_path = __file__[:-11] + 'lib' + platform.architecture()[0][:2] + '/v' + sys.version[:3] + '/'
if lib_path not in sys.path:
    sys.path.append(lib_path)
print(sys.path)

from awGA import awGA
from awGA_templet import awGA_templet
from bs2int import bs2int
from bs2rv import bs2rv
from crtbase import crtbase
from crtbp import crtbp
from crtfld import crtfld
from crtip import crtip
from crtpp import crtpp
from crtrp import crtrp
from frontplot import frontplot
from i_awGA_templet import i_awGA_templet
from indexing import indexing
from migrate import migrate
from mut import mut
from mutate import mutate
from mutbga import mutbga
from mutbin import mutbin
from mutgau import mutgau
from mutint import mutint
from mutpp import mutpp
from ndomindeb import ndomindeb
from ndomin import ndomin
from ndominfast import ndominfast
from nsga2_templet import nsga2_templet
from powing import powing
from q_sorted_templet import q_sorted_templet
from ranking import ranking
from recdis import recdis
from recint import recint
from reclin import reclin
from recombin import recombin
from reins import reins
from rep import rep
from rwGA import rwGA
from rws import rws
from scaling import scaling
from selecting import selecting
from sga_code_templet import sga_code_templet
from sga_new_code_templet import sga_new_code_templet
from sga_permut_templet import sga_permut_templet
from sga_new_permut_templet import sga_new_permut_templet
from sga_real_templet import sga_real_templet
from sga_new_real_templet import sga_new_real_templet
from sus import sus
from tour import tour
from trcplot import trcplot
from upNDSet import upNDSet
from xovdp import xovdp
from xovdprs import xovdprs
from xovmp import xovmp
from xovpm import xovpm
from xovsh import xovsh
from xovshrs import xovshrs
from xovsp import xovsp
from xovsprs import xovsprs
