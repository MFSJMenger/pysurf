# -*- coding: utf-8 -*-

"""Console script for pysurf."""
import sys
import click


@click.command("Hallo Welt")
@click.option('--db', default='db', help='this is a database')
def main(db=4):
    """Console script for pysurf."""
    click.echo("Replace this message by putting your code into "
               "pysurf.cli.main")
    click.echo(type(db))
    click.echo("See click documentation at https://click.palletsprojects.com/")
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
