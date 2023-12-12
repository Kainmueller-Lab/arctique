import pip
import os

from pathlib import Path

cur_dir = str(Path(os.path.dirname(os.path.realpath(__file__))).parent.parent.resolve())

pip.main(['install', cur_dir, '--user'])