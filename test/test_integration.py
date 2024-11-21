import os
from pathlib import Path

cur_path = Path(os.path.dirname(os.path.realpath(__file__)))


class TestIntegrationSuite:

    def test_main(self):
        from render import main
        import sys
        # set args
        sys.argv = [
            "render.py",
            "--start-idx=42",
            "--n-samples=1"
        ]

        # run main
        main()