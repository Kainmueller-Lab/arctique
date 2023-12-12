import argparse

import rHEnder


def startup():
    """Entry points of `polarityjam`."""
    parser = RHEnderParser().parser
    args = parser.parse_args()

    __run_entrypoint(args, parser)


def __run_entrypoint(args, parser):
    """Calls a specific subcommand."""
    pass


class ArgumentParser(argparse.ArgumentParser):
    """Override default error method of all parsers to show help of sub-command."""

    def error(self, message: str):
        self.print_help()
        self.exit(2, '%s: error: %s\n' % (self.prog, message))


class RHEnderParser(ArgumentParser):
    def __init__(self):
        super().__init__()
        self.parser = self.create_parser()

    @staticmethod
    def create_parser() -> ArgumentParser:
        """Parent parser for all subparsers to have the same set of arguments."""
        parser = ArgumentParser(add_help=True)
        parser.add_argument('--version', '-V', action='version', version="%s " % rHEnder.__version__)
        # parse logging
        parser.add_argument(
            '--log',
            required=False,
            help='Logging level for your command. Choose between %s' %
                 ", ".join(['INFO', 'DEBUG']),
            default='INFO'
        )
        parser.add_argument(
            '--myarg1', '-a', required=False, help='My first argument', default="Hello World"
        )

        return parser
