import logging
import argparse
from argtools import Command


class RecipeBook(Command):
    def __init__(self):
        self.env = Env()
        super(RecipeBook, self).__init__()

    def __call__(self, cmd, name=None):
        """
        Use this as decorator.

        # example
        book = RecipeBook()

        @book
        def cmd1(g):
            g.do('something')
        """
        return self.add_recipe(cmd, name=name)

    def add(self, cmd):
        return self.call(cmd)

    def add_recipe(self, cmd, name=None):
        name = cmd.__name__ if name is None else name
        new = Command()
        new._take_over(cmd)
        if name:
            new.__name__ = name

        self._children.append(new)
        recipe = Recipe(new, parent=self)
        return recipe

    def _add_parser_rules(self, parser):
        # TODO exclusive ?
        parser.add_argument('--run', action='store_true')
        parser.add_argument('--show', action='store_true')

        super(RecipeBook, self)._add_parser_rules(parser)

    def run(self, args=None):
        parser = argparse.ArgumentParser(description=self.description,
                                         epilog=self.epilog,
                                         formatter_class=self.formatter_class)
        # setup root parser
        self._add_parser_rules(parser)
        print((self._children))
        print((vars(self.env)))
        # setup subcommand parsers
        if self._children:
            subs = parser.add_subparsers(help='valid subcommands')
            for subcmd in self._children:
                #subcmd = recipe.get_command()
                if subcmd.name is None:
                    continue # skip trailing subcommand
                subp = subs.add_parser(subcmd.name,
                                       description=subcmd.description,
                                       epilog=subcmd.epilog,
                                       help=subcmd.description,
                                       formatter_class=self.formatter_class)
                subcmd._add_parser_rules(subp)

        args = parser.parse_args(args=args)
        for before_run in self._before_runs:
            before_run(args)

        has_func = hasattr(args, 'func')
        if not has_func:
            parser.error('too few arguments')
            sys.exit(1)

        return self._run(args.func, args)

    def _run(self, builder, args):
        self._setup_logger(args)

        graph = TaskGraph(env=self.env)
        for key, val in list(vars(args).items()):
            setattr(graph.env, key, val)

        builder(graph)
        if args.run:
            print(('mode is run', graph))
        else:
            print(('mode is show', graph))

        # shortname = os.path.basename(sys.argv[0])
        # try:
        #     pass
        # except IOError as e:
        #     if e.errno != 32:  # ignore SIGPIPE
        #         raise
        # finally:
        # self.logger.info('start %s', shortname)
        # self.logger.info('end %s', shortname)


class Recipe:
    def __init__(self, cmd, parent=None):
        assert callable(cmd)
        if not isinstance(cmd, Command):  # wrap by Command
            cmd = Command()(cmd)
        self._cmd = cmd
        self._parent = parent

    def get_command(self):
        return self._cmd

    def __call__(self, *args, **kwds):   # build recipe here
        parser = argparse.ArgumentParser(description=self._cmd.description,
                                 epilog=self._cmd.epilog,
                                 formatter_class=self._cmd.formatter_class)
        self._cmd._add_parser_rules(parser)
        parser.set_defaults(**kwds)
        args = parser.parse_args(args=args)
        graph = TaskGraph(env=self._parent.env)
        args.func(graph)
        # build recipe
        return graph



class Env(object):
    def __init__(self, parent=None):
        super(Env, self).__init__()
        self._parent = parent

    def __getattr__(self, key):
        try:
            return getattr(self, key)
        except KeyError:
            return getattr(self._parent, key)


class TaskGraph:
    def __init__(self, env=None):
        self._tasks = []
        self.env = Env(parent=env)
        self._check_list = []
        self._make_list = []
        self._info_list = []

    def do(self, *args, **kwds):
        self._tasks.append(('do', args, kwds))  # TODO

    def check(self, fstr):
        self._check_list.append(fstr)

    def make(self, fstr):
        self._make_list.append(fstr)

    def info(self, fstr):
        self._info_list.append(fstr)

    def set(self, key, value, format=True):
        self.env.key = (value if not format else value.format(self.env))


# decorator
def argument(*args, **kwds):
    return Command().add_argument(*args, **kwds)

#recipe = RecipeBook()
