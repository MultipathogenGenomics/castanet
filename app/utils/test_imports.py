
from app.utils.error_handlers import stoperr


def import_test():
    try:
        import numpy
        import pandas
        import matplotlib
        import pydantic
        import fastapi
        import uvicorn
        import autopep8
        import reportlab
        import biotite
        import plotly
        import rapidfuzz
        import Bio
        import kaleido
    except Exception as ex:
        stoperr(
            f"Couldn't find Python libraray: {ex}. Ensure you're running the correct version of Python and libraries are installed (pip install -r requirements.txt)")
