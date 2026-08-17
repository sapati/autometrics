from pathlib import Path
import typer
from .config import load_config
from .pipeline import run as run_pipeline

app = typer.Typer(help="Automated variable-star CCD photometry")


@app.command()
def run(project: Path = typer.Argument(..., exists=True, file_okay=False), target: str = typer.Argument(...), config: Path | None = typer.Option(None, "--config", "-c")):
    """Calibrate, solve, identify, measure, and report one observing run."""
    try:
        run_pipeline(load_config(project, config), target)
    except Exception as exc:
        typer.echo(f"autometrics: {exc}", err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
