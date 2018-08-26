""" It's a test of importing geatpy"""

import sys
import platform

lib_path = __file__[:-11] + 'lib' + platform.architecture()[0][:2] + '/v' + sys.version[:3] + '/'
if lib_path not in sys.path:
    sys.path.append(lib_path)
print(sys.path)
print()
import geatpy
