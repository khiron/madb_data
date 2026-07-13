import time
from pathlib import Path

figures_path = Path("..") / "figures"

class pdf_writer:
    """class that handles super annooying mathjax warning box in plotly pdf's"""

    def __init__(self) -> None:
        self._done_once = False
        figures_path.mkdir(parents=True, exist_ok=True)

    def __call__(self, fig, filename: str) -> None:
        # the sleep, plus successive write, is ESSENTIAL to avoid the super annoying
        # "[MathJax]/extensions/MathMenu.js" text box error
        # but we only need to do this once
        
        # Ensure filename has .pdf suffix
        path = Path(filename)
        if path.suffix.lower() != ".pdf":
            path = path.with_suffix(".pdf")

        full_path = figures_path / path
        if not self._done_once:
            fig.write_image(full_path)
            time.sleep(2)
            self._done_once = True
        fig.write_image(full_path)