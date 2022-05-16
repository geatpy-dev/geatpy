import argparse

import geatpy

parser = argparse.ArgumentParser()
parser.add_argument('-v',
                    '--version',
                    help='show version',
                    action='version',
                    version='geatpy {}'.format(geatpy.__version__))
parser.parse_args()
