import os
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

