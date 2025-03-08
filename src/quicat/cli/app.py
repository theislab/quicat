import click
from rich.console import Console
from rich.traceback import install

from quicat.cli.extract import Extractor
from quicat.cli.simulation import Simulator

# Install rich traceback for better error formatting
install()
console = Console()


@click.group()
def cli() -> None:
    """quicat: A CLI tool for barcode processing and analysis."""
    pass


@click.command()
@click.argument("config_path", type=click.Path(exists=True))
def extract(config_path: str) -> None:
    """
    Extract barcodes using the configuration file provided in CONFIG_PATH.
    """
    try:
        extractor = Extractor(config_path)
        extractor.run()
        console.print("[bold green]Barcode extraction completed successfully![/bold green]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise


@click.command()
@click.argument("config_path", type=click.Path(exists=True))
def simulate(config_path: str) -> None:
    """
    Simulate barcodes using the configuration file provided in CONFIG_PATH.
    """
    try:
        simulator = Simulator(config_path)
        simulator.generate_synthetic_data()
        console.print("[bold green]Simulation completed successfully![/bold green]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise


cli.add_command(extract)
cli.add_command(simulate)

if __name__ == "__main__":
    cli()
