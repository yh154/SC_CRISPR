from __future__ import absolute_import

# Import happens on setup, but on setup networks dependencies won't have
# been installed. So if we can't import, don't worry about it.]
try:
    from counts import UMIClusterer
except ImportError:
    pass