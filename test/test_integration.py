import os
from pathlib import Path

cur_path = Path(os.path.dirname(os.path.realpath(__file__)))


class TestIntegrationSuite:

    def test_main(self):
        import main  # runs the main script
