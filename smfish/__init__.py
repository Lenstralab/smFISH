from . import FISH_pipeline, spotDetection_functions
from ._version import __version__, __git_commit_hash__

# Adding local tllab_common to path if it exists
import sys
import os
tlpath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'tllab_common')
if os.path.exists(tlpath):
    sys.path.insert(0, tlpath)
